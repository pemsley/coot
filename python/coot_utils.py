# coot-utils.py
# adapted from coot-utils.scm
#
# Copyright 2004, 2005, 2006, 2007 by Bernhard Lohkamp
# Copyright 2008, 2009 by Bernhard Lohkamp, The University of York
# Copyright 2000 by Paul Emsley
# Copyright 2004, 2005, 2006, 2007 by Paul Emsley, The University of York
#    <one line to give the program's name and a brief idea of what it does.>
#    Copyright (C) <year>  <name of author>
#    Copyright 2013, 2014, 2016 by Medical Research Council
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

import os
import re
import string
import numbers
import math
import urllib

import coot
# import coot_gui # circular dependency
from coot_redefine_functions import * # 20220828-PE this needs fixing

# hack this in for now
global use_gui_qm
use_gui_qm = False


# 3D annotations - a bit of a hack currently
global annotations
annotations = []

# thank you ebassi!
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


# used in Extensions -> Representation -> Ball & Stick
global default_ball_and_stick_selection
default_ball_and_stick_selection = "//A/1-2"

# for mingw debug
global have_mingw
have_mingw = False
if os.getenv("MSYSTEM"):
    have_mingw = True

# this is set by the main application: False or 1 (use this python gui)
# or 2 (guile-gtk is being used so use this if gui-gtk function is
# not available)

global user_defined_alert_smarts
# example: user_defined_alert_smarts = [['C', 'my-user-defined alert for carbon']]
user_defined_alert_smarts = []


# not sure if the implementation of the macros will work

# 'Macro' to tidy up a a setup of functions to be run with no backup
# for a particular molecule.
#
# funcs is a normal set of functions (not a thunk), here i.e. a list of
# functions with function as a list with func name and args,
# e.g.: [centre_of_mass, 0], [func_name, arg1, arg2,...],...
#
def with_no_backups(imol, *funcs):
    """Execute functions with backups disabled for a molecule.
    
    Temporarily disables backups for the specified molecule, runs the provided
    functions, then restores the original backup state.
    
    Parameters
    ----------
    imol : int
        Molecule number
    *funcs : tuple
        Variable number of function calls, each as [func, arg1, arg2, ...]
        
    Examples
    --------
    >>> with_no_backups(0, [coot.refine_zone, 0, "A", 10, 15, ""])
    """

    b_state = coot.backup_state(imol)
    coot.turn_off_backup(imol)
    for f in funcs:
        func = f[0]
        args = f[1:len(f)]
        # print "BL DEBUG:: func %s and args %s" %(func, args)
        func(*args)
    if backup_mode == 1:
        coot.turn_on_backup(imol)


# 'Macro' to tidy up a set of functions to be run with automatic
# accepting of the refinement
# returns the result of last function run...
#
# funcs is a normal set of functions (not a thunk), here i.e. a list of
# functions with function as a list with func name and args,
# e.g.: [centre_of_mass, 0], [func_name, arg1, arg2,...],...
#
def with_auto_accept(*funcs):
    """Execute functions with automatic acceptance of refinement.
    
    Temporarily enables immediate replacement mode, runs the provided functions
    with automatic acceptance, then restores the original refinement state.
    
    Parameters
    ----------
    *funcs : tuple
        Variable number of function calls, each as [func, arg1, arg2, ...]
        
    Returns
    -------
    object
        Result of the last function executed
        
    Examples
    --------
    >>> with_auto_accept([coot.refine_zone, 0, "A", 10, 15, ""])
    """

    replace_state = coot.refinement_immediate_replacement_state()
    coot.set_refinement_immediate_replacement(1)
    ret = None
    for f in funcs:
        func = f[0]
        args = f[1:len(f)]
        try:
            func(*args)
            ret = coot.accept_moving_atoms_py()
        except SystemError as e:
            print("Py -----------------------------")
            print('WARNING:: with_auto_accept() caught a SystemError exception')
            print(e)
            print("Py -----------------------------")

    if (replace_state == 0):
        coot.set_refinement_immediate_replacement(0)

    return ret   # returns result of last functions!!!!

# 'Macro' to run funcs on an active atom
# funcs is function, active_atom specifiers and extra args
# func, args, "aa_imol", "aa_chain_id", ..., args
# or list thereof
# [[func1, extra_arg1, ..., "aa_imol", "aa_chain_id",..., extra_arg2, extra arg3, ..], [func2,...]]
# returns what? The value from the last function evaluated
#


def using_active_atom(*funcs):
    """Execute function using the currently active (closest) atom.
    
    Finds the closest atom to the screen center and passes its specifications
    to the provided function. Use placeholder strings like "aa_imol", "aa_chain_id",
    etc. to reference atom properties.
    
    Parameters
    ----------
    *funcs : tuple
        Function and arguments with placeholders for active atom specs.
        Available placeholders:
        - "aa_imol": molecule number
        - "aa_chain_id": chain ID
        - "aa_res_no": residue number
        - "aa_ins_code": insertion code
        - "aa_atom_name": atom name
        - "aa_alt_conf": alternate conformation
        - "aa_res_spec": residue spec [chain_id, res_no, ins_code]
        
    Examples
    --------
    >>> using_active_atom(coot.refine_zone, "aa_imol", "aa_chain_id", 
    ...                   "aa_res_no", "aa_res_no", "")
    """

    active_atom = closest_atom_simple()
    print("debug:: active_atom", active_atom)
    if not active_atom:
        coot.add_status_bar_text("No residue found")
    else:
        def convert_arg(arg):
            try:
                c_arg = aa_dict[arg]
                return c_arg
            except KeyError as e:
                return arg

        aa_dict = {"aa_imol":      active_atom[0],
                   "aa_chain_id":  active_atom[1],
                   "aa_res_no":    active_atom[2],
                   "aa_ins_code":  active_atom[3],
                   "aa_atom_name": active_atom[4],
                   "aa_alt_conf":  active_atom[5],
                   "aa_res_spec":  [active_atom[1],  # chain_id
                                    active_atom[2],  # res_no
                                    active_atom[3]]} # ins_code

        if True:
            func = funcs[0]
            args = funcs[1:]
            c_args = [convert_arg(item) for item in args]
            # print("debug:: func", func)
            # print("debug:: c_args", c_args)
            func(*c_args)


# #  This is what using_active_atom() used to be - I don't like it
            # if isinstance(item, list):
            #     arg_ls = []
            #     for ele in item:
            #         arg_ls.append(arg_to_append(ele))
            #     return arg_ls
            # else:
            #     if item in aa_dict:
            #         return aa_dict[item]
            #     else:
            #         return item

#         if (len(funcs) == 1):
#             # we have a list of functions
#             # so use
#             ls_funcs = funcs[0]
#         elif (type(funcs[0]) is ListType):
#             # we have a range of lists
#             # use as is
#             ls_funcs = funcs
#         else:
#             # we have a single function with args
#             # make into list
#             ls_funcs = [funcs]

#         for ele in ls_funcs:
#             func = ele[0]
#             func_args = ele[1:]
#             args = []
#             for arg in func_args:
#                 ins = arg_to_append(arg)
#                 args.append(ins)
#             ret = func(*args)
#         return ret

# here some truly pythonic version of the macros. Should replace
# them in usage too:

# Pythonic 'Macro' to tidy up a a setup of functions to be run with no backup
# for a particular molecule.
#
# use with 'with', e.g.:
#
# > with NoBackups(imol=0):
#      coot.refine_zone(imol, "A", 43, 45, "")
#      coot.accept_regularizement()
#


class NoBackups:
    """'Macro' to tidy up a a setup of functions to be run with no backup
    for a particular molecule (default imol=0).

    use with 'with', e.g.:

    > with WithNoBackups(imol=0):
        coot.refine_zone(imol, "A", 43, 45, "")
        coot.accept_regularizement()
    """

    def __init__(self, imol=0):
        self.imol = imol

    def __enter__(self):
        self.b_state = coot.backup_state(self.imol)
        coot.turn_off_backup(self.imol)

    def __exit__(self, type, value, traceback):
        if (self.b_state == 1):
            coot.turn_on_backup(self.imol)

# Pythonic 'Macro' to tidy up a set of functions to be run with automatic
# accepting of the refinement.
#
#    use with 'with', e.g.:
#
#    >with AutoAccept():
#        coot.refine_zone(0, "A", 43, 45, "")
#


class AutoAccept:
    """
    Pythonic 'Macro' to tidy up a set of functions to be run with automatic
    accepting of the refinement.

    use with 'with', e.g.:

    > with AutoAccept():
         coot.refine_zone(0, "A", 43, 45, "")

    """

    def __init__(self):
        self.replace_state = -1
        pass

    def __enter__(self):
        self.replace_state = coot.refinement_immediate_replacement_state()
        coot.set_refinement_immediate_replacement(1)

    def __exit__(self, type, value, traceback):
        coot.accept_regularizement()
        if (self.replace_state == 0):
            coot.set_refinement_immediate_replacement(0)


class UsingActiveAtom:
    """
    Run functions on the active atom.

    use with 'with', e.g.:

    > with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf]:
          coot.refine_zone(aa_imol, aa_chain_id, aa_res_no-2, aa_res_no+2, aa_ins_code)

    alternative usage to get res_spec as well

    > with coot_utils.UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf, aa_res_spec]:
          coot.refine_zone(aa_imol, aa_chain_id, aa_res_no-2, aa_res_no+2, aa_ins_code)

    """

    def __init__(self, with_res_spec=False):
        self.no_residue = False
        self.res_spec = with_res_spec
        pass

    def __enter__(self):
        self.active_atom = active_residue()
        if (not self.active_atom):
            coot.add_status_bar_text("No (active) residue found")
            self.no_residue = True
            #self.__exit__(None, "dummy", None)
        else:
            imol      = self.active_atom[0]
            chain_id  = self.active_atom[1]
            res_no    = self.active_atom[2]
            ins_code  = self.active_atom[3]
            atom_name = self.active_atom[4]
            alt_conf  = self.active_atom[5]
            res_spec  = [self.active_atom[1],  # chain_id
                         self.active_atom[2],  # res_no
                         self.active_atom[3]] # ins_code
            if self.res_spec:
                return [imol, chain_id, res_no, ins_code, atom_name, alt_conf, res_spec]
            else:
                return [imol, chain_id, res_no, ins_code, atom_name, alt_conf]

    def __exit__(self, type, value, traceback):
        if (self.no_residue):
            # internal calling of exit, ignore errors
            return True
        pass


def coot_home():
    """Return the normalized HOME directory (or COOT_HOME on Windows).
    If not available, return False"""
    import os
    global have_mingw

    home = os.getenv('HOME')
    if (os.name == 'nt' and not have_mingw):
        home = os.getenv('COOT_HOME')
    if not home:
        # no home!? Poor you.
        return False
    else:
        return os.path.normpath(home)


# make the directory and return the directory name. If you cant make
# it directly in this directory try to make it in $HOME. Return False
# if complete failure. e.g. coot-ccp4 or coot-backup
#
def get_directory(dir_name):
    """make the directory and return the directory name. If you cant make
    it directly in this directory try to make it in $HOME. Return False
    if complete failure. e.g. coot-ccp4 or coot-backup"""

    import os
    if os.path.isdir(dir_name):
        return dir_name
    if os.path.isfile(dir_name):
        return False
    else:
        status = coot.make_directory_maybe(dir_name)
        if (status == 0):
            return dir_name
        else:
            h = coot_home()
            if not h:
                # failed to get home!? Mmmh
                return False
            else:
                new_dir = os.path.join(h, dir_name)
                status = coot.make_directory_maybe(new_dir)
                if (status == 0):
                    return new_dir
                else:
                    return False

# 20220228-PE accepted on merge - maybe wrong.
def get_active_molecule():
    """Get the first displayed model molecule number.
    
    Returns
    -------
    int or False
        Molecule number of first displayed model, or False if none found
    """
    active_atom = closest_atom_simple()
    if active_atom:
        return active_atom[0]
    else:
        return -1


# Pythonize function: return a python boolean.
#
def molecule_has_hydrogens(imol):
    """Check if molecule contains hydrogen atoms.
    
    Parameters
    ----------
    imol : int
        Molecule number
        
    Returns
    -------
    bool
        True if molecule contains hydrogens, False otherwise
    """
    return (coot.molecule_has_hydrogens_raw(imol) == 1)


def add_hydrogens_using_refmac(imol):
    """Add hydrogen atoms to molecule using REFMAC.
    
    Parameters
    ----------
    imol : int
        Molecule number to add hydrogens to
    """
    out_file_name = os.path.join("coot-refmac", molecule_name_stub(imol, 0) + '-needs-H.pdb')
    in_file_name  = os.path.join("coot-refmac", molecule_name_stub(imol, 0) +  '-with-H.pdb')
    coot.make_directory_maybe('coot-refmac')
    coot.write_pdb_file(imol, out_file_name)
    return add_hydrogens_using_refmac_inner(imol, in_file_name, out_file_name)


def add_hydrogens_to_chain_using_refmac(imol, chain_id):
    out_file_name = os.path.join("coot-refmac",
                                 molecule_name_stub(imol, 0) + '-chain-' + chain_id + '-needs-H.pdb')
    in_file_name = os.path.join("coot-refmac",
                                molecule_name_stub(imol, 0) + '-chain-' + chain_id + '-with-H.pdb')
    coot.make_directory_maybe('coot-refmac')
    coot.write_chain_to_pdb_file(imol, chain_id, out_file_name)
    return add_hydrogens_using_refmac_inner(imol, in_file_name, out_file_name)


def add_hydrogens_using_refmac_inner(imol, in_file_name, out_file_name):

    status = popen_command("refmac5",
                           ['XYZIN', out_file_name, 'XYZOUT', in_file_name],
                           ['MAKE HOUT YES', 'NCYCLE 0', 'END'],
                           'refmac-H-addition.log', 0)
    try:
        if status == 0:
            # all good
            return coot.add_hydrogens_from_file(imol, in_file_name)
        else:
            return False
    except:
        return False


# set this to a function accepting two argument (the molecule number
# and the manipulation mode) and it will be run after a model
# manipulation.
# e.g.
# def post_manipulation_script(imol, mode):
#  ... continue below
#
# The manipulation mode will be one of (MOVINGATOMS), (DELETED) or
# (MUTATED) and these can be tested with "=".
#
# e.g.
#
# if (mode == DELETED):
#      display/print "Something was deleted"
#
# This is a global variable now, so that it can be used within other functions
#
global post_manipulation_script
post_manipulation_script = False

# similar for the active residue
# do something based on the active residue (presumably)
#
global post_set_rotation_centre_script
post_set_rotation_centre_script = False

global post_read_model_hook
post_read_model_hook = False

# return a boolean
#


def pre_release_qm():
    return "-pre" in coot.coot_version()


def molecule_list(molecule_filter_function):
    ls = []
    for mol_no in range(coot.graphics_n_molecules()):
        # print("######## debug molecule_list()", mol_no, molecule_filter_function(mol_no))
        if molecule_filter_function(mol_no):
            ls.append(mol_no)
    # print("######### debug molecule_list() returning ls", ls)
    return ls

# return a list of molecule numbers (closed and open)
# The elements of the returned list need to be tested against
# is_valid_model_molecule_qm
#
def molecule_number_list():
    """Get list of all valid molecule numbers.
    
    Returns
    -------
    list of int
        List of all valid molecule numbers (both models and maps)
    """
    ret = []
    for mol_no in range(coot.graphics_n_molecules()):
        if valid_map_molecule_qm(mol_no) or valid_model_molecule_qm(mol_no):
            ret.append(mol_no)
    return ret


def model_molecule_number_list():
    """Get list of all model molecule numbers.
    
    Returns
    -------
    list of int
        List of all valid model (not map) molecule numbers
    """
    return list(filter(valid_model_molecule_qm, molecule_number_list()))

def display_all_maps():
    map_list = map_molecule_list()
    # 20220420-PE this map doesn't work - I don't understand why
    # map(lambda imol: coot.set_map_displayed(imol, 1), map_list)
    for m in map_list:
        coot.set_map_displayed(m, 1)

# c.f. graphics_info_t::undisplay_all_model_molecules_except(int imol)


def undisplay_all_maps_except(imol_map):

    map_list = map_molecule_list()
    for imol in map_list:
        if (imol != imol_map):
            coot.set_map_displayed(imol, 0)
    coot.set_map_displayed(imol_map, 1)

# 20220228-PE was called "just_one_or_next" I think
def display_cycle_through_maps():

    def next_map(current_map_number, map_number_list):
        try:
            current_idx = map_number_list.index(current_map_number)
        except:
            current_idx = -1
        l = len(map_number_list)
        if current_idx > -1:
            next_index = 0 if current_idx + 1 == l else current_idx + 1
            return map_number_list[next_index]
        return map_number_list[0]

    # main line

    map_list = map_molecule_list()
    current_displayed_maps = [imol for imol in map_list if coot.map_is_displayed(imol) == 1]
    n_displayed = len(current_displayed_maps)

    if not current_displayed_maps:
       if map_list:
          undisplay_all_maps_except(map_list[0])
    else:
       if n_displayed == 1:
          cm = current_displayed_maps[0]
          if cm == map_list[-1]:
             # show all
             display_all_maps()
          else:
             n = next_map(cm, map_list)
             undisplay_all_maps_except(n)
       else:
          undisplay_all_maps_except(map_list[0])

# is isn't quite because one of the options is "all"
# legacy for installed key bindings:
def just_one_or_next_map():
    return display_cycle_through_maps()


# Test for prefix-dir (1) being a string (2) existing (3) being a
# directory (4-6) modifiable by user (ie. u+rwx).  prefix_dir must be a
# string.
#
# Return True or False.
#
def directory_is_modifiable_qm(prefix_dir):
    from types import StringType
    ret = False
    # check string:
    ret = type(prefix_dir) is StringType
    if ret:
        # check existence:
        ret = os.access(prefix_dir, os.F_OK)
        if ret:
            # check dir:
            ret = os.path.isdir(prefix_dir)
            if ret:
                # check readability
                ret = os.access(prefix_dir, os.R_OK)
                if ret:
                    # check writability
                    ret = os.access(prefix_dir, os.W_OK)
                    if ret:
                        # check executability (needed?!)
                        ret = os.access(prefix_dir, os.X_OK)
    return ret

# return an absolute file-name for file-name or False
#
# def absolutify(file_name)  - exists in python as os.path.abspath use that!!

# Find the most recently created file from the given glob and dir
#
# return False on no-such-file
#


def most_recently_created_file(glob_str, dir):

    import glob
    import time

    patt = os.path.join(dir, glob_str)
    files = glob.glob(patt)

    latest_file = False
    latest_mtime = 0

    for file_ in files:
        this_mtime = os.path.getmtime(file_)
        if this_mtime > latest_mtime:
            latest_file = file_
            latest_mtime = this_mtime

    return latest_file


# Convert a residue_spec to an mmdb atom selection string.
# FIXME:: to be tested
#
def residue_spec_to_atom_selection_string(centre_residue_spec):
    ret = "//" + centre_residue_spec[0] + \
          "/" + str(centre_residue_spec[1])
    return ret


def residue_atom_to_atom_name(ra):
    if not isinstance(ra, list):
        return False
    else:
        return ra[0][0]

def residue_atom_to_position(ra):
    if not isinstance(ra, list):
        return False
    else:
        return ra[2]


#         I don't like these function names - remove them at some stage

def residue_atom2atom_name(ra):
    if not isinstance(ra, list):
        return False
    else:
        return ra[0][0]


def residue_atom2alt_conf(ra):
    if not isinstance(ra, list):
        return False
    else:
        return ra[0][1]


def residue_atom2occupancy(ra):
    if not isinstance(ra, list):
        return False
    else:
        return ra[1][0]


def residue_atom2position(ra):
    if not isinstance(ra, list):
        return False
    else:
        return ra[2]

# start using new convention, maybe. could do some tests if atom name and
# alt conf is not False
# residue_info atom needs other parameters to make a spec for an atom


def residue_atom_to_atom_spec(ra, chain_id, res_no, ins_code):
    if not isinstance(ra, list):
        return False
    else:
        return [chain_id, res_no, ins_code,
                residue_atom2atom_name(ra),
                residue_atom2alt_conf(ra)]


def residue_spec_to_chain_id(rs):
    """Extract chain ID from residue spec.
    
    Parameters
    ----------
    rs : list
        Residue spec [chain_id, res_no, ins_code]
        
    Returns
    -------
    str
        Chain ID
        
    Examples
    --------
    >>> spec = ["A", 42, ""]
    >>> residue_spec_to_chain_id(spec)
    'A'
    """
    if not isinstance(rs, list):
        return False
    else:
        if (len(rs) == 3):
            return rs[0]
        else:
            if (len(rs) == 4):
                return rs[1]
            else:
                return False


def residue_spec_to_res_no(rs):
    """Extract residue number from residue spec.
    
    Parameters
    ----------
    rs : list
        Residue spec [chain_id, res_no, ins_code]
        
    Returns
    -------
    int
        Residue number
        
    Examples
    --------
    >>> spec = ["A", 42, ""]
    >>> residue_spec_to_res_no(spec)
    42
    """
    if not isinstance(rs, list):
        return False
    else:
        if (len(rs) == 3):
            return rs[1]
        else:
            if (len(rs) == 4):
                return rs[2]
        return False


def residue_spec_to_ins_code(rs):
    """Extract insertion code from residue spec.
    
    Parameters
    ----------
    rs : list
        Residue spec [chain_id, res_no, ins_code]
        
    Returns
    -------
    str
        Insertion code
        
    Examples
    --------
    >>> spec = ["A", 42, "A"]
    >>> residue_spec_to_ins_code(spec)
    'A'
    """
    if not isinstance(rs, list):
        return False
    else:
        if (len(rs) == 3):
            return rs[2]
        else:
            if (len(rs) == 4):
                return rs[3]
        return False


def residue_specs_match_qm(spec_1, spec_2):
    """Check if two residue specs refer to the same residue.
    
    Parameters
    ----------
    spec_1 : list
        First residue spec [chain_id, res_no, ins_code]
    spec_2 : list
        Second residue spec [chain_id, res_no, ins_code]
        
    Returns
    -------
    bool
        True if specs match, False otherwise
    """
    if (residue_spec_to_chain_id(spec_1) ==
            residue_spec_to_chain_id(spec_2)):
        if (residue_spec_to_res_no(spec_1) ==
                residue_spec_to_res_no(spec_2)):
            if (residue_spec_to_ins_code(spec_1) ==
                    residue_spec_to_ins_code(spec_2)):
                return True
    return False


def atom_spec_to_imol(atom_spec):
    """Extract molecule number from atom spec.
    
    Parameters
    ----------
    atom_spec : list
        Atom spec [imol, chain_id, res_no, ins_code, atom_name, alt_conf, ...]
        
    Returns
    -------
    int
        Molecule number
    """
    import types
    if not (isinstance(atom_spec, list)):
        return False
    else:
        if (len(atom_spec) == 6):
            return atom_spec[0]
        if (len(atom_spec) == 7):
            return atom_spec[1]
        return False


def residue_spec_to_residue_name(imol, spec):
    if not isinstance(spec, list):
        return False
    if (len(spec) == 3):
        return coot.residue_name(imol,
                            spec[0],
                            spec[1],
                            spec[2])
    elif (len(spec) == 4):
        return coot.residue_name(imol,
                            spec[1],
                            spec[2],
                            spec[3])
    else:
        return False

# for sorting residue specs
#


def residue_spec_less_than(spec_1, spec_2):
    chain_id_1 = residue_spec_to_chain_id(spec_1)
    chain_id_2 = residue_spec_to_chain_id(spec_2)
    if chain_id_2 < chain_id_1:
        return True
    else:
        rn_1 = residue_spec_to_res_no(spec_1)
        rn_2 = residue_spec_to_res_no(spec_2)
        if rn_2 < rn_1:
            return True
        else:
            ins_code_1 = residue_spec_to_ins_code(spec_1)
            ins_code_2 = residue_spec_to_ins_code(spec_2)
            return ins_code_2 < ins_code_1


def residue_spec_to_string(spec):
    """Convert residue spec to human-readable string.
    
    Parameters
    ----------
    spec : list
        Residue spec [chain_id, res_no, ins_code]
        
    Returns
    -------
    str
        Formatted string like "A42" or "B15A" (with insertion code)
        
    Examples
    --------
    >>> residue_spec_to_string(["A", 42, ""])
    'A42'
    >>> residue_spec_to_string(["B", 15, "A"])
    'B15A'
    """
    ret = residue_spec_to_chain_id(spec) + " "
    ret += str(residue_spec_to_res_no(spec))
    ret += residue_spec_to_ins_code(spec)
    return ret

# Return a list of molecules that are maps
#


def map_molecule_list():
    """Get list of all map molecule numbers.
    
    Returns
    -------
    list of int
        List of all valid map molecule numbers
        
    Examples
    --------
    >>> maps = map_molecule_list()
    >>> print(f"Found {len(maps)} maps")
    """

    map_list = []
    for i in range(coot.graphics_n_molecules()):
        if coot.is_valid_map_molecule(i) == 1:
            map_list.append(i)
    return map_list

# Return a list of molecules that are (coordinate) models
#


def model_molecule_list():
    """Get list of all model molecule numbers (same as model_molecule_number_list).
    
    Returns
    -------
    list of int
        List of all valid model molecule numbers
    """

    model_list = []
    for i in range(coot.graphics_n_molecules()):
        if coot.is_valid_model_molecule(i) == 1:
            model_list.append(i)
    return model_list

# Return True(False) if @var{imol} is (isn't) a shelx molecule.
#


def shelx_molecule_qm(imol):
    """Check if molecule is a SHELX molecule.
    
    Parameters
    ----------
    imol : int
        Molecule number
        
    Returns
    -------
    bool
        True if molecule is SHELX format, False otherwise
    """
    if (coot.is_shelx_molecule(imol) == 1):
        return True
    else:
        return False


# return an int. 0 means no, 1 means yes, -1 on error
#
is_protein_chain_qm = coot.is_protein_chain_p

# Is a nucleotide chain?
# Now return a boolean
#


def is_nucleotide_chain_qm(imol, chain_id):
    """Check if chain is nucleotide/nucleic acid.
    
    Parameters
    ----------
    imol : int
        Molecule number
    chain_id : str
        Chain ID
        
    Returns
    -------
    bool
        True if chain is nucleotide, False otherwise
    """
    return coot.is_nucleotide_chain_p(imol, chain_id) == 1


# Set the virtual trackball behaviour.
#
# trackball @var{type} is a string: either 'flat' or 'spherical-surface'
#
def set_virtual_trackball_type(type):
    if (type == "flat"):
        coot.vt_surface(1)
    elif (type == "spherical-surface"):
        coot.vt_surface(0)
    else:
        print("virtual trackball type", type, "not understood")

# Is ls a list of strings? Return True or False
#


def list_of_strings_qm(ls):
    import types
    not_str = 0
    if type(ls) is not ListType:
        return False
    else:
        for item in ls:
            if isinstance(item, (str,)):
                pass
            else:
                not_str += 1
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
    """Get current rotation centre coordinates.
    
    Returns
    -------
    list of float
        [x, y, z] coordinates of rotation centre
    """
    return [coot.rotation_centre_position(0),
            coot.rotation_centre_position(1),
            coot.rotation_centre_position(2)]

# this is actually not essentail since python has these funtion(s)


def number_list(a, b):
    result = []
    if a == b:
        result.append(a)
        return result
    elif a > b:
        return result
    else:
        while a <= b:
            result.append(a)
            a = a + 1
        return result


def file_n_lines(file_name):
    if not os.path.isfile(file_name):
        return False
    else:
        fin = open(file_name, 'r')
        n_lines = sum(1 for line in fin)
        fin.close()
        return n_lines

# backport, so that we can replace once move to Python 2.7 is done


def check_output(*popenargs, **kwargs):
    r"""Run command with arguments and return its output as a byte string.

Backported from Python 2.7 as it's implemented as pure python on stdlib.

>>> check_output(['/usr/bin/python', '--version'])
Python 2.6.2
"""

    import subprocess
    process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        error = subprocess.CalledProcessError(retcode, cmd)
        error.output = output
        raise error
    return output

# returns empty string if there is a problem running cmd
#


def shell_command_to_string(cmd):

    import subprocess
    import sys

    major, minor, micro, releaselevel, serial = sys.version_info

    if (major >= 2 and minor >= 7):
        try:
            ret = subprocess.check_output(cmd.split())
        except:
            ret = ""
    else:
        try:
            ret = check_output(cmd)
        except:
            ret = ""
    return ret

# Return True or False
# adapted from find_exe
# this finds absolute file names too
#


def command_in_path_qm(cmd, no_disk_search=True,
                       only_extension_here=None,
                       add_extensions_here=[]):

    """Check if a command is available in the system PATH.

    Parameters
    ----------
    cmd : str
        Command name to search for
    no_disk_search : bool, optional
        If True, only check PATH without disk search (default: True)

    Returns
    -------
    str or False
        Full path to command if found, False otherwise

    Examples
    --------
    >>> command_in_path_qm("refmac5")
    '/usr/local/bin/refmac5'
    """

    exe_path = find_exe(cmd, "PATH", no_disk_search=True,
                        screen_info=False,
                        only_extension=only_extension_here,
                        add_extensions=add_extensions_here)

    if exe_path:
        return True
    else:
        return False


def command_in_path_qm_old_version(cmd, only_extension="", add_extensions=[]):
    # test for command (see goosh-command-with-file-input description)
    #
    import os
    import string

    # we shall check for full path names first
    if (os.path.isfile(cmd)):
        return True
    else:
        extensions = []
        cmd_noext = cmd
        # some windows magic
        if (os.name == 'nt'):
            file_ext = file_name_extension(cmd)
            cmd_noext = strip_extension(cmd)
            if (file_ext):
                extensions = [file_ext]
            else:
                tmp_ext = os.environ["PATHEXT"].split(os.pathsep)
                # list of extensions (no dot) only
                extensions = [ext[1:] for ext in tmp_ext]

        if only_extension:
            extensions = [only_extension]
        if add_extensions:
            extensions += add_extensions

        program_names = [cmd_noext]
        if extensions:
            program_names += [cmd_noext + "." + ext for ext in extensions]

        try:
            primary_path = os.environ["PATH"]
            for cmd_name in program_names:
                for path in string.split(primary_path, os.pathsep):
                    program_exe = os.path.join(path, cmd_name)
                    #           print "BL DEBUG:: program_exe is", program_exe
                    if (os.path.isfile(program_exe)):
                        return True
            return False
        except:
            print("BL WARNING:: couldnt open $PATH")  # this shouldnt happen
            return False


global gtk_thread_return_value
gtk_thread_return_value = None
# Where cmd is e.g. "refmac"
#       args is ["HKLIN","thing.mtz"]
#       data_list is ["HEAD","END"]
#       log_file_name is "refmac.log"
#       screen_flag True/False to display or not in shell window
#       stderr_capture True/False pipe stderr to stdout too
#       local_env can be set to change the environment variables the
#                 command is run in.
#
# Return the exist status e.g. 0 or 1. Or False if cmd not found.
#
# uses os.popen if python version < 2.4 otherwise subprocess
#


def popen_command(cmd, args, data_list, log_file, screen_flag=False,
                  stderr_capture=False, local_env=None):

    import sys
    import string
    import os
    import subprocess

    major, minor, micro, releaselevel, serial = sys.version_info

    if os.path.isfile(cmd):
        cmd_execfile = cmd
    else:
        if not(command_in_path_qm(cmd)):
            print("command ", cmd, " not found in $PATH!")
            print("BL INFO:: Maybe we'll find it somewhere else later...")
        cmd_execfile = find_exe(cmd, "CBIN", "CBIN", "CCP4_BIN", "PATH")

    if cmd_execfile:

        if True:

            print("debug open write", log_file)

            log = open(log_file, 'wb')

            cmd_args = [cmd_execfile] + args

            # set stdin
            pipe_data=False
            if isinstance(data_list, str) and os.path.isfile(data_list):
                stdin_inp = open(data_list)
            else:
                pipe_data=True
                stdin_inp = subprocess.PIPE

            # set stderr
            if stderr_capture:
                stderr_arg = subprocess.STDOUT
            else:
                stderr_arg = None

            if screen_flag:
                process = subprocess.Popen(cmd_args, stdin=stdin_inp,
                                           stdout=subprocess.PIPE,
                                           stderr=stderr_arg,
                                           env=local_env)
            else:
                process = subprocess.Popen(cmd_args, stdin=stdin_inp,
                                           stdout=log, stderr=stderr_arg,
                                           env=local_env)

            if pipe_data:
                for data in data_list:
                    # process.stdin.write(data + "\n")
                    process.stdin.write('{}\n'.format(data).encode('utf-8'))
                process.stdin.close()
            else:
                # file?
                if isinstance(stdin_inp, file):
                    stdin_inp.close()

            if screen_flag:
                for line in process.stdout:
                    # remove trailing whitespace
                    # print("#", line.rstrip(" \n"))
                    print("#", line)
                    log.write(line)
            process.wait()
            log.close()

            return process.returncode

        else:
            # popen (old)
            data_list_file = "data_list_file_tmp.txt"

            # make args string
            args_string = string.join(args)

            # write tmp input file
            input = file(data_list_file, 'w')
            for data in data_list:
                input.write(data + '\n')
            input.close()

            # print "BL DEBUG:: popen command is", cmd_execfile + " " + args_string + " < " + data_list_file + " > " + log_file
            status = os.popen(cmd_execfile + " " + args_string +
                              " < " + data_list_file + " > " + log_file, 'r')
            cmd_finished = status.close()
            if (cmd_finished != None):
                print("running command ", cmd, " failed!")
                return 1
            else:
                # print "cmd ", cmd ," seems to have run ok!"
                # log_file_size = os.stat(log_file)[stat.ST_SIZE]
                return 0
            os.remove(data_list_file)
    else:
        print("WARNING:: could not find %s, so not running" % cmd)
        return False


def ok_popen_status_qm(status):

    if not isNumber(status):
        return False
    else:
        return status == 0


# example usage:
# popen_command("mtzdump",["HKLIN","a.mtz"],["HEAD","END"],"test.log",0)

# Crude test to see of 2 floats are the same (more or less).
# Used in a unit test after setting the atom position.
#
def close_float_qm(x1, x2):

    return (abs(x1 - x2) < 0.001)

# "a.b.res" -> "a.b"
# file_name_sans_extension
#


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

    import os
    import string
    root, ext = os.path.splitext(file_name)
    if ext:
        ext = string.lstrip(ext, '.')
        return ext
    else:
        return ""

# e.g. "a.pdb" -> "a-tmp.pdb"
#


def add_tmp_extension_to(file_name):

    import types
    if isinstance(file_name, (str,)):
        root, ext = os.path.splitext(file_name)
        f = root + "-tmp" + ext
        return f
    else:
        return "tmp"

# Same function as strip_extension, different name, as per scsh, in fact.
#


def file_name_sans_extension(s):
    return strip_extension(s)

# /a/b.t -> b.t   d/e.ext -> e.ext
# file-name-sans-path


def strip_path(s):
    import os
    head, tail = os.path.split(s)
    return tail

# does s start with a "/" ?
# return True or False
#
# for windows return True when drive letter, e.g. C, or even \\ (backslash):


def slash_start_qm(s):
    import types
    import string
    if isinstance(s, (str,)):
        if len(s) > 0:
            if (s.startswith("/") or s.startswith("\\")):
                return True
            else:
                if (os.name == 'nt' and len(s) > 1):
                    # for windows check if the start is a drive letter
                    drive_letter_ls = [
                        dl + ':' for dl in string.ascii_uppercase]
                    if (s[0:2] in drive_letter_ls):
                        return True
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
#


def every_nth(ls, n):

    elements = list(range(0, len(ls), n))
    a = []
    for i in elements:
        a.append(ls[i])
    return a

# now usefull in testing
#
# residue_atoms must be a list
#


def get_atom_from_residue(atom_name, residue_atoms, alt_conf):
    """Get atom_info from a residue.
    residue_atoms must be a list
    """

    if ((isinstance(residue_atoms, list)) and
            (residue_atoms != [])):
        for residue_atom in residue_atoms:
            if (residue_atom[0][0] == atom_name and
                    residue_atom[0][1] == alt_conf):
                return residue_atom
    # print "BL WARNING:: no atom name %s found in residue" %atom_name
    return False  # no residue name found fail save

#


def get_atom_from_spec(imol, atom_spec):
    return get_atom(imol,
                    atom_spec_to_chain_id(atom_spec),
                    atom_spec_to_res_no(atom_spec),
                    atom_spec_to_ins_code(atom_spec),
                    atom_spec_to_atom_name(atom_spec),
                    atom_spec_to_alt_loc(atom_spec))

# return atom info or False (if atom not found).
#


def get_atom(imol, chain_id, resno, ins_code, atom_name, alt_conf_internal=""):

    res_info = coot.residue_info_py(imol, chain_id, resno, "")

    if (not res_info):
        return False
    else:
        ret = get_atom_from_residue(atom_name, res_info, alt_conf_internal)
        return ret

#


def residue_info_dialog_displayed_qm():
    if (coot.residue_info_dialog_is_displayed() == 1):
        return True
    else:
        return False

# multi_read pdb reads all the files matching
# @code{@emph{glob_pattern}} in
# directory @code{@emph{dir}}.  Typical usage of this might be:
# @code{multi_read_pdb("a*.pdb",".")}
# BL says: in windows dir needs the 'C:/' pattern, '/c/'won't work
#


def multi_read_pdb(glob_pattern, dir):
    import glob
    import os
    patt = os.path.normpath(dir+'/*.'+glob_pattern)
    all_files = glob.glob(patt)
    for file in all_files:
        print("BL INFO:: reading ", file)
        coot.read_pdb(file)

# read_pdb_all reads all the "*.pdb" files in the current directory.
#


def read_pdb_all():
    import glob
    import os
    recentre_status = coot.recentre_on_read_pdb()
    coot.set_recentre_on_read_pdb(0)
    patt = os.path.normpath(os.path.abspath(".")+'/*.pdb')
    all_files = glob.glob(patt)
    for file in all_files:
        print("BL INFO:: reading ", file)
        coot.read_pdb(file)
    coot.set_recentre_on_read_pdb(recentre_status)

# return the dir-name on success.
#
# return False if dir_name is a file or we can't do the mkdir.
#


def coot_mkdir(dir_name):
    """return the dir-name on success.

    return False if dir_name is a file or we can't do the mkdir."""
    import os

    if (os.path.isfile(dir_name)):
        return False
    else:
        if (os.path.isdir(dir_name)):
            return dir_name
        else:
            try:
                os.mkdir(dir_name)
                return dir_name
            except:
                # failed to do a mkdir
                # better to return false
                return False

# return the view matrix (useful for molscript, perhaps).
# BL says: like all matrices is a python list [...]
#


def view_matrix():
    return [get_view_matrix_element(row_number, column_number) for row_number in range(3) for column_number in range(3)]

# return the transposed view matrix (useful for molscript, perhaps).
# BL says: like all matrices is a python list [...]
#


def view_matrix_transp():
    return [get_view_matrix_element(column_number, row_number) for row_number in range(3) for column_number in range(3)]

# return the view quaternion
#


def view_quaternion():

    ret = list(map(coot.get_view_quaternion_internal, [0, 1, 2, 3]))
    return ret

# Return the view number
#


def add_view(position, quaternion, zoom, view_name):

    args = position + quaternion
    args.append(zoom)
    args.append(view_name)
    ret = coot.add_view_raw(*args)
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
        elif (x > 0 and y < 0):
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

    ls = list(map(convert_sign, pr[1:], [m21 - m12, m02 - m20, m10 - m01]))

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

    coot.set_view_quaternion(matrix2quaternion(m00, m10, m20,
                                          m01, m11, m21,
                                          m02, m12, m22))

# Miguel's molecular orientation axes
#


def miguels_axes():
    coot.set_axis_orientation_matrix(*view_matrix())
    coot.set_axis_orientation_matrix_usage(1)

# Return the molecule centre as a list of 3 numbers.
#
#  Note: mol_cen could contain values less than -9999.
#


def molecule_centre(imol):
    return [coot.molecule_centre_internal(imol, 0),
            coot.molecule_centre_internal(imol, 1),
            coot.molecule_centre_internal(imol, 2)]

# Move the centre of molecule number imol to the current screen centre
#


def move_molecule_to_screen_centre(imol):
    if valid_model_molecule_qm(imol):
        rotate_centre = rotation_centre()
        coot.translate_molecule_by(imol, (rotate_centre[0]-molecule_centre(imol)[0]),
                              (rotate_centre[1]-molecule_centre(imol)[1]),
                              (rotate_centre[2]-molecule_centre(imol)[2]))


# This is a short name for the above.
# deftexi move_molecule_here
move_molecule_here = move_molecule_to_screen_centre

# Return a nine-membered list of numbers.
#


def identity_matrix():
    return [1, 0, 0, 0, 1, 0, 0, 0, 1]

# e.g. translation('x',2)
#  -> [2, 0, 0]
# Return: False on error
#


def translation(axis, length):
    import operator
# BL says: we dont check if axis is string, yet at least not directly
    if (isinstance(length, numbers.Number)):
        if (axis == "x"):
            return [length, 0, 0]
        elif (axis == "y"):
            return [0, length, 0]
        elif (axis == "z"):
            return [0, 0, length]
        else:
            print("symbol axis: ", axis, " incomprehensible")
            return False
    else:
        print("incomprehensible length argument: ", length)
        return False

# Rotate degrees about screen axis, where axis is either 'x', 'y' or 'z'.
#


def rotate_about_screen_axis(axis, degrees):
    import math
    import operator

    def deg_to_rad(degs):
        return (degs * 3.1415926 / 180.0)

    def simple_rotation_x(alpha):
        cos_alpha = math.cos(alpha)
        sin_alpha = math.sin(alpha)
        return [1, 0, 0, 0, cos_alpha, -sin_alpha, 0, sin_alpha, cos_alpha]

    def simple_rotation_y(alpha):
        cos_alpha = math.cos(alpha)
        sin_alpha = math.sin(alpha)
        return [cos_alpha, 0, sin_alpha, 0, 1, 0, -sin_alpha, 0, cos_alpha]

    def simple_rotation_z(alpha):
        cos_alpha = math.cos(alpha)
        sin_alpha = math.sin(alpha)
        return [cos_alpha, -sin_alpha, 0, sin_alpha, cos_alpha, 0, 0, 0, ]

# BL says: I dont know what the next 2 defines are for...
# looks not used and/or useless to me
# seems that only 2nd matrix is used and not view_matrix!
    def vm(): view_matrix()
    def mult(mat1, mat2): mat2
# end of uselessness...

    if (isinstance(degrees, numbers.Number)):
        if (axis == "x"):
            mult(view_matrix(), simple_rotation_x(deg_to_rad(degrees)))
        elif (axis == "y"):
            mult(view_matrix(), simple_rotation_y(deg_to_rad(degrees)))
        elif (axis == "z"):
            mult(view_matrix(), simple_rotation_z(deg_to_rad(degrees)))
        else:
            print("symbol axis: ", axis, " incomprehensible")
    else:
        print("incomprehensible length argument: ", degrees)


# Support for old toggle functions.  (consider instead the raw
# functions use the direct set_displayed functions).
#
def toggle_display_map(imol, idummy):
    if (coot.map_is_displayed(imol) == 0):
        coot.set_map_displayed(imol, 1)
    else:
        coot.set_map_displayed(imol, 0)

# toggle the display of imol
#


def toggle_display_mol(imol):
    if (coot.mol_is_displayed(imol) == 0):
        coot.set_mol_displayed(imol, 1)
    else:
        coot.set_mol_displayed(imol, 0)

# toggle the active state (clickability) of imol
#


def toggle_active_mol(imol):
    if (coot.mol_is_active(imol) == 0):
        coot.set_mol_active(imol, 1)
    else:
        coot.set_mol_active(imol, 0)

# return a python (list) representation of molecule imol, or False if we can't
# do it (imol is a map, say)
# optional arg: chain
#


def python_representation(imol, chains=[]):

    if (not valid_model_molecule_qm(imol)):
        return False
    else:
        ls = []

        def r_info(imol, chain_id, n):
            res_name = coot.resname_from_serial_number(imol, chain_id, n)
            res_no = coot.seqnum_from_serial_number(imol, chain_id, n)
            ins_code = coot.insertion_code_from_serial_number(imol, chain_id, n)
            return [res_no, ins_code, res_name, coot.residue_info_py(imol, chain_id, res_no, ins_code)]

        if not chains:
            # use all
            chains = chain_ids(imol)
        ls = [[[chain_id, [r_info(imol, chain_id, serial_number) for serial_number in range(
            coot.chain_n_residues(chain_id, imol))]] for chain_id in chains]]
        return ls

# reorder chains
#


def reorder_chains(imol):

    # reorder elements of chain_list: e.g.
    #
    # chain_list: [["C", [xx]], ["A", [xx]], ["B", [xx]]]
    #
    def reorder_chains_in_model(chain_list):
        list(map(lambda model: model.sort(), chain_list))

    p_rep = python_representation(imol)
    if (type(p_rep) is ListType):
        reorder_chains_in_model(p_rep)
        clear_and_update_molecule(imol, p_rep)


# transform a coordinates molecule by a coot-rtop (which is a Python
# expression of a clipper::RTop), i.e. a list of a 9-element list and
# a 3 element list, e.g. [[1, 0, 0, 0, 1, 0, 0, 0, 1], [4.5, 0.4, 1.2]]
#
def transform_coords_molecule(imol, rtop):

    ls = []
    for i in rtop:
        for j in i:
            ls.append(j)

    coot.transform_molecule_by(imol, *ls)

# @code{transform_map(imol, mat, trans, about_pt, radius, space_group, cell)}
#
# where space_group is a HM-symbol and cell is a list of 6
# parameters, where the cell angles are in degrees.
#
# or @code{transform_map(imol, trans, about_pt, radius)} for a simple translation
#
# or @code{transform_map(imol, trans, radius)} when using the default
# rotation-centre as the about-pt
#
# returns new map mol number or None if no map could be transformed/created
#


def transform_map(*args):

    ret = None

    def tf(imol, mat, trans, about_pt, radius, space_group, cell):

        print("DEBUG:: here in tf with ", imol, mat, trans, about_pt, radius, space_group, cell)

        return transform_map_raw(imol,
                                 mat[0], mat[1], mat[2],
                                 mat[3], mat[4], mat[5],
                                 mat[6], mat[7], mat[8],
                                 trans[0], trans[1], trans[2],
                                 about_pt[0], about_pt[1], about_pt[2],
                                 radius,
                                 space_group,
                                 cell[0], cell[1], cell[2],
                                 cell[3], cell[4], cell[5])

    print("len args is", len(args))

    # main line
    if (len(args) == 7):
        ret = tf(args[0], args[1], args[2], args[3], args[4], args[5], args[6])
    # imol_map mat trans about_pt radius:
    elif (len(args) == 5):
        imol = args[0]
        ret = tf(imol, args[1], args[2], args[3], args[4],
                 space_group(imol), cell(imol))
    # no matrix specified:
    elif len(args) == 4:
        imol = args[0]
        print("calling tf with ", imol, identity_matrix(), args[1], args[2], args[3], space_group(imol), cell(imol))
        ret = tf(imol, identity_matrix(), [args[1], args[2], args[3]],
                 rotation_centre(), cell(imol)[0],
                 space_group(imol), cell(imol))
    # no matrix or about point specified:
    elif (len(args) == 3):
        imol = args[0]
        ret = tf(args[0], identity_matrix(), args[1], rotation_centre(),
                 args[2], space_group(imol), cell(imol))
    else:
        print("arguments to transform-map incomprehensible: args: ", args)
    return ret


# return then NCS master of the first molecule that has ncs.
#
# return "" on fail to find an ncs chain
#
def get_first_ncs_master_chain():

    r = ""
    for mols in model_molecule_list():
        ncs_masters = ncs_master_chains(mols)
        if ncs_masters:
            return ncs_masters[0]
    return r

# Define a map transformation function that obeys Lapthorn's Law of
# NCS Handling Programs
#
# typical usage: transform_map_using_lsq_matrix(1, "A", 10, 30, 0, "A", 10, 30, 2, rotation_centre(), 6)
#
# Remember, that now the about-pt is the "to" point, i.e. the maps are brought from
# somewhere else and generated about the about-pt.
#


def transform_map_using_lsq_matrix(imol_ref, ref_chain, ref_resno_start, ref_resno_end,
                                   imol_mov, mov_chain, mov_resno_start, mov_resno_end,
                                   imol_map, about_pt, radius):

    coot.clear_lsq_matches()
    coot.add_lsq_match(ref_resno_start, ref_resno_end, ref_chain,
                  mov_resno_start, mov_resno_end, mov_chain, 1)
    space_group = symmetry_operators_to_xHM(symmetry_operators(imol_ref))
    cell_params = cell(imol_ref)

    if not (space_group and cell):
        message = "Bad cell or symmetry for molecule" + \
            str(cell) + str(space_group) + str(imol_ref)
        print("Bad cell or symmetry for molecule", message)
        ret = -1           # invalid mol! or return message!?

    else:
        rtop = apply_lsq_matches(imol_ref, imol_mov)
        ret = transform_map(imol_map, rtop[0], rtop[1], about_pt, radius,
                            space_group, cell_params)
    return ret

# Make the imol-th map brighter.
#
# Scale_factor > 1 makes brighter...
#

# brighten maps
def brighten_map(imol, scale_factor):

    if valid_map_molecule_qm(imol):
        current_colour = coot.map_colour_components_py(imol)
        try:
            new_v = []
            for i in range(len(current_colour)):
                new_v.append(current_colour[i] * float(scale_factor))
                if new_v[i] < 0.05:
                    new_v[i] = 0.05
                elif new_v[i] > 1.0:
                    new_v[i] = 1.0
                else:
                    pass
            coot.set_map_colour(imol, *new_v)
        except:
            print("bad non-list current-colour ", current_colour)
        coot.graphics_draw()

# Make all maps brighter
#
def brighten_maps():
    list(map(lambda imap: brighten_map(imap, 1.25), map_molecule_list()))

# Make all maps darker
#
def darken_maps():
    list(map(lambda imap: brighten_map(imap, 0.8), map_molecule_list()))

# pastelize the colour of the map by degree
def pastelize_map(imol_map: int, degree: float):
    if valid_map_molecule_qm(imol_map):
        colour = coot.map_colour_components_py(imol_map)
        try:
            for i in range(len(colour)):
                cc = colour[i];
                r = 1.0 - cc # extra "whiteness" (for this colour component)
                colour[i] += r * degree;
                colour[i] *= (1-0 - 0.5 * degree) # I don't want bright pastels
            coot.set_map_colour(imol_map, *colour)
        except Exception as e:
            print("ERROR:: pastelize_map(): ", e, "bad non-list current-colour ", colour)
        coot.graphics_draw()

# pastelize maps
def pastelize_maps():
    list(map(lambda imap: pastelize_map(imap, 0.25), map_molecule_list()))

def solidify_maps():
    list(map(lambda imap: coot.set_draw_solid_density_surface(imap, 1), map_molecule_list()))

def lineify_maps():
    list(map(lambda imap: coot.set_draw_map_standard_lines(imap, 1), map_molecule_list()))


# return a list of chain ids for given molecule number @var{imol}.
# return empty list on error
#
def chain_ids(imol):

    """ Get list of chain IDs for a molecule

    Parameters
    ----------
    imol : int
        molecule number

    Returns
    -------
    list
        list of chain ID strings (e.g., ['A', 'B', 'C'])
    """

    # if this fails, check that you have not set chain_id to override coot's chain_id()
    return [chain_id(imol, ic) for ic in range(coot.n_chains(imol))]

# convert from interface name to schemisch name to be equivalent to Paul's naming
#
# return True or False
#


def is_solvent_chain_qm(imol, chain_id):
    if coot.is_solvent_chain_p(imol, chain_id) == 1:
        return True
    else:
        return False

# convert from interface name to schemisch name to be equivalent to Paul's naming
#
# return True or False
#


def is_protein_chain_qm(imol, chain_id):
    return coot.is_protein_chain_p(imol, chain_id) == 1

# convert from interface name to schemisch name to be equivalent to Paul's naming
#
# return True or False
#


def is_nucleotide_chain_qm(imol, chain_id):
    """Check if chain is nucleotide/nucleic acid.
    
    Parameters
    ----------
    imol : int
        Molecule number
    chain_id : str
        Chain ID
        
    Returns
    -------
    bool
        True if chain is nucleotide, False otherwise
    """
    return coot.is_nucleotide_chain_p(imol, chain_id) == 1

# python (schemeyish) interface to eponymous scripting interface function!?
# return True or False
#


def valid_model_molecule_qm(imol):
    """Check if molecule number is a valid model.
    
    Parameters
    ----------
    imol : int
        Molecule number to check
        
    Returns
    -------
    bool
        True if imol is a valid model molecule, False otherwise
    """
    if (coot.is_valid_model_molecule(imol) == 1):
        return True
    else:
        return False

# python (schemeyish) interface to eponymous scripting interface function.
# return True or False
#


def valid_map_molecule_qm(imol):
    """Check if molecule number is a valid map.
    
    Parameters
    ----------
    imol : int
        Molecule number to check
        
    Returns
    -------
    bool
        True if imol is a valid map molecule, False otherwise
    """
    if coot.is_valid_map_molecule(imol) == 1:
        return True
    else:
        return False

# python (schemeyish) interface to eponymous scripting interface function.
# return True or False
#
def valid_map_with_associated_data_molecule_qm(imol):
    if coot.is_valid_map_molecule(imol) == 1:
        p = coot.refmac_parameters_py(imol)
        if p:
            return True
        else:
            return False
    else:
        return False

# convenience function (slightly less typing).
#
# Return True or False
#


def valid_refinement_map_qm():
    return valid_map_molecule_qm(coot.imol_refinement_map())

# python (schemeyish) interface to shelx molecule test
#
# Return True or False
#


def shelx_molecule_qm(imol):
    """Check if molecule is a SHELX molecule.
    
    Parameters
    ----------
    imol : int
        Molecule number
        
    Returns
    -------
    bool
        True if molecule is SHELX format, False otherwise
    """
    return coot.is_shelx_molecule(imol) == 1

# python (schemeyish) interface to the function that returns whether or not a map
# is a difference map.
#
# Return True or False.
#


def is_difference_map_qm(imol_map):
    if not valid_map_molecule_qm(imol_map):
        return False
    else:
        return coot.map_is_difference_map(imol_map) == 1

# Does residue resno with insertion code ins_code of chain chain_id
# and in molecule number imol exist?
#
# Return True or False
#


def residue_exists_qm(imol, chain_id, resno, ins_code):
    return coot.does_residue_exist_p(imol, chain_id, resno, ins_code) == 1

# Does the residue contain hetatoms?
# Return True or False.
#


def residue_has_hetatms_qm(imol, chain_id, res_no, ins_code):
    return coot.residue_has_hetatms(imol, chain_id, res_no, ins_code) == 1

# Return a list of 3 float for the centre of mas of molecule number imol.
#
# on failure return False.
#


def centre_of_mass(imol):
    centre = eval(coot.centre_of_mass_string(imol))
    if (centre == 0):
        print("molecule number", imol, "is not valid")
        return False
    else:
        print("Centre of mass for molecule %s is %s" % (imol, centre))
        # for use somewhere let's return the centre
        return centre

# Return as a list the occupancy, temperature_factor, x y z coordinates
# of the given atom.
# (the x,y,z are in Cartesian Angstroms).
#
# on error (e.g. atom not found) return False
#


def atom_specs(imol, chain_id, resno, ins_code, atom_name, alt_conf):

    return atom_info_string(imol, chain_id, resno, ins_code, atom_name, alt_conf)


def atom_spec_to_string(spec):

    ret = atom_spec_to_chain_id(spec) + " "
    ret += str(atom_spec_to_res_no(spec))
    ret += atom_spec_to_ins_code(spec) + " "
    ret += atom_spec_to_atom_name(spec)
    al = atom_spec_to_alt_loc(spec)
    if al:
        ret += " " + al
    return ret


def atom_spec_to_residue_spec(atom_spec):
    l = len(atom_spec)
    if l == 5:
        return atom_spec[:3]
    else:
        if l == 6:  # active_residue give an atom-spec prepended by the imol
            return atom_spec[1:][:3]
        else:
            if l == 7:  # active_residue give an atom-spec prepended a bool and by the imol. These should die
                spec = atom_spec[2:][:3]
                print("atom_spec to residue_spec returns", spec)
                return spec
            return None

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
            if coot.map_is_difference_map(map_mol) == 0:
                return map_mol
        print("BL WARNING:: we couldnt find a non difference map for fitting!")
        return -1


# Ian Tickle says (as far as I can understand) that the target rmsd
# should be 0.25 or thereabouts.  You can over-ride it now.
#
# BL says:: still, I dont get it to 0.25, so we make it 1.0
# Even setting it to 0.8 (more or less arbitratry, but accounting for
# some redundancy in the parameters (Ian's comment)) is not enough, I think
# It seems reality is different to theory, here at least...
global target_auto_weighting_value
target_auto_weighting_value = 1.0

# Set the refinement weight (matrix) by iterating the refinement and
# varying the weight until the chi squares (not including the
# non-bonded terms) reach 1.0 =/- 10%.  It uses sphere refinement.
# The refinement map must be set!!  At the end show the new weight in
# the status bar.  Seems to take about 5 rounds. (bails out after 20)
#


def auto_weight_for_refinement():

    global target_auto_weighting_value
    # return a pair of the imol and a list of residue specs.
    # or False if that is not possible

    def sphere_residues(radius):
        active_atom = active_residue()
        if not active_atom:    # check for list?
            print("No active atom")
            return False
        else:
            centred_residue = active_atom[1:4]
            imol = active_atom[0]
            other_residues = coot.residues_near_residue_py(imol, centred_residue, radius)
            all_residues = [centred_residue]
            if (isinstance(other_residues, list)):
                all_residues += other_residues
            return [imol, all_residues]

    # the refinement function that is run and returns nice refinement
    # results
    #
    def refinement_func():
        sr = sphere_residues(3.5)
        if sr:
            ret = with_auto_accept([refine_residues, sr[0], sr[1]])
            return ret
            # return with_auto_accept([refine_residues, sr[0], sr[1]])
        else:
            return False

    # get rid of non-bonded chi-squared results from the input list ls.
    #
    def no_non_bonded(ls):
        ret = []
        if ls:
            for item in ls:
                if not (item[0] == "Non-bonded"):
                    ret.append(item)
        return ret

    # return False or a number, which is the current overweighting of the
    # density terms.  (of course, when not overweighted, the geometric
    # chi squareds will be about 1.0).
    #
    def weight_scale_from_refinement_results(rr):
        if not rr:   # check for list?
            return False
        else:
            results_inner = rr[2]
            if not results_inner:
                return False
            else:
                nnb_list = no_non_bonded(results_inner)
                chi_squares = [x[2] for x in nnb_list]
                n = len(chi_squares)
                summ = sum(chi_squares)
                if n == 0:
                    return False
                else:
                    return summ/n

    # main body
    #
    results = refinement_func()
    n_trials = 0
    while results:
        av_rms_d = weight_scale_from_refinement_results(results)
        print("av_rms_d:", av_rms_d)
        if not av_rms_d:   # check for number?
            return False
        if n_trials > 20:
            print("BL INFO:: refinement did not converge, bailing out")
            return False
        if (av_rms_d < (target_auto_weighting_value * 1.1) and
                av_rms_d > (target_auto_weighting_value * 0.9)):
            # done
            s = "Success: Set weight matrix to " + str(coot.matrix_state())
            coot.add_status_bar_text(s)
            break
        else:
            # more refinement required

            # squared causes ringing,
            # as does 1.5.
            # Simple is overdamped.
            current_weight = coot.matrix_state()
            new_weight = (target_auto_weighting_value *
                          current_weight) / av_rms_d
            print("INFO:: setting refinement weight to %s from * %s / %s"
                  % (new_weight, current_weight, av_rms_d))
            if (new_weight < 2):
                # weight refinement not converging
                print("BL INFO:: not convering, weight to set was", new_weight)
            coot.set_matrix(new_weight)
            results = refinement_func()
        n_trials += 1

    return True   # return successful termination...

# Print the sequence of molecule number @var{imol}
#
# This is not really a util, perhaps it should be somewhere else?
#


def print_sequence(imol):

    for chain in chain_ids(imol):
        coot.print_sequence_chain(imol, chain)

# simple utility function to return the contents of a file as a string.
#


def pir_file_name2pir_sequence(pir_file_name):
    import os
    if (not os.path.isfile(pir_file_name)):
        return False
    else:
        try:
            fin = open(pir_file_name, 'r')
            str = fin.read()
            fin.close()
            return str
        except:
            return False

# Associate the contents of a PIR file with a molecule.
#


def associate_pir_file(imol, chain_id, pir_file_name):
    seq_text = pir_file_name2pir_sequence(pir_file_name)
    if seq_text:
        coot.assign_pir_sequence(imol, chain_id, seq_text)
    else:
        print("WARNING:: associate-pir-file: bad text for", pir_file_name)

# Associate the contents of a fasta file with a molecule.
#


def associate_fasta_file(imol, chain_id, pir_file_name):
    seq_text = pir_file_name2pir_sequence(pir_file_name)
    if seq_text:
        coot.assign_fasta_sequence(imol, chain_id, seq_text)
    else:
        print("WARNING:: associate-fasta-file: bad text for", pir_file_name)


# comma key hook
def graphics_comma_key_pressed_hook():
    pass

# dot key hook


def graphics_dot_key_pressed_hook():
    pass

# a list of [code, key, name, thunk]
# e.g. [103, "g", "Goto Blob", coot.blob_under_pointer_to_screen_centre()]


global key_bindings
# we shall see if it exists, if not initialize it
try:
    key_bindings
except:
    key_bindings = []


def decode_key(key_val_name):
    try:
        import gtk
        key_value = int(gtk.gdk.keyval_from_name(key_val_name))
        # on some windows: special characters seem to have high value,
        # so need to convert these properly too
        if (not key_value or key_value >= 100000):
            # new python needs a long there!? I hope it wont harm old!?
            new_val = int(ord(key_val_name))
            key_value = int(gtk.gdk.unicode_to_keyval(new_val))
        return key_value
    except:
        return key_sym_code(key_val_name)

# Add a key binding
#
# with a given name, key (e.g. "x" or "S") and the function to run
# (a thunk) when that key is pressed.
#


def add_key_binding(name, key, thunk):

    if use_gui_qm:

        from types import IntType, StringType

        global key_bindings, std_key_bindings
        std_keys = [elem[1] for elem in std_key_bindings]
        keys = [elem[1] for elem in key_bindings]
        codes = [elem[0] for elem in key_bindings]
        if (key in std_keys):
            print("INFO:: you shall not overwrite a standard key binding (%s)" % key)
        else:
            if (type(key) is IntType):
                if (key in keys):
                    print("INFO:: you are overwriting existing key", key)
                    key_bindings.pop(keys.index(key))
                key_bindings.append([key, key, name, thunk])
            elif (type(key) is StringType):
                code = decode_key(key)
                if (code in codes):
                    print("INFO:: you are overwriting existing key (from code)", key)
                    key_bindings.pop(codes.index(code))
                if "Control_" in key:
                    code = decode_key(key[-1])
                if (("Control_" in key) or not (code == -1)):
                    key_bindings.append([code, key, name, thunk])
                else:
                    print("INFO:: key %s not found in code table" % key)
            else:
                print("BL WARNING:: invalid key", key)


# general key press hook, not for public use!!
#
def graphics_general_key_press_hook(key, control_flag=0):
    global key_bindings
    # print "graphics_general_key_press_hook(): Key %s was pressed" %key
    if control_flag:
        codes = [elem[0] for elem in key_bindings if "Control_" in elem[1]]
        funcs = [elem[3] for elem in key_bindings if "Control_" in elem[1]]
    else:
        codes = [elem[0] for elem in key_bindings if not "Control_" in elem[1]]
        funcs = [elem[3] for elem in key_bindings if not "Control_" in elem[1]]
    if (key in codes):
        index = codes.index(key)
        func = funcs[index]
        # print "BL DEBUG:: index and executing:", index, func
        apply(func)
    else:
        if coot.coot_has_guile() and is_windows():
            coot.run_scheme_command("(graphics-general-key-press-hook " +
                               str(key) +
                               ")")
        print("Key %s not found in (python) key bindings" % key)


# Function requested by Mark White.
#
# read XtalView (and maybe other) .vu files and convert them into generic
# objects.
#
# Pass the filename and an object name e.g.
# read_vu_file("axes.vu", "axes")
#
# Returns: nothing interesting.
#
def read_vu_file(filename, obj_name):

    from types import StringType
    from operator import isNumberType

    def colour_from_number(obj):
        if type(obj) is StringType:
            return obj
        elif isNumbertype(obj):
            if (obj == 0):
                return "white"
            else:
                return "red"
        else:
            return "white"

    # main body
    n = coot.new_generic_object_number(obj_name)
    fin = open(filename, 'r')
    lines = fin.readlines()
    fin.close()
    for line in lines:
        current_line = line.split()
        if len(current_line) == 7:
            colour = colour_from_number(current_line[-1])
            try:
                coords = list(map(float, current_line[0:-1]))
            except:
                print("BL WARNING:: cannot make float from cordinates",
                      current_line[0:-1])
                return
            coot.to_generic_object_add_line(n, colour, 2,
                                       *coords)
    coot.set_display_generic_object(n, 1)


# residue_test_func is a function that takes 4 arguments, the
# chain_id, resno, inscode and residue_serial_number (should it be
# needed) and returns either False or return something interesting
# (e.g. text for a button label).
#
# Return a list of residues, each of which has a return value at the
# start, ie. [return_value, chain_id, res_no, ins_code]
#
def residues_matching_criteria(imol, residue_test_func):

    matchers = []
    # these specs are prefixed by the serial number
    for molecule_residue_specs in all_residues_with_serial_numbers(imol):
        rs = molecule_residue_specs[1:]
        if residue_test_func(residue_spec_to_chain_id(rs),
                             residue_spec_to_res_no(rs),
                             residue_spec_to_ins_code(rs),
                             molecule_residue_specs[0]):
            matchers.append(rs)
    return matchers

# Now this is in the API
#
# Return residue specs for all residues in imol (each spec is preceeded by True)
#


def all_residues(imol):
    r = all_residues_with_serial_numbers(imol)
    try:
        return [e[1:] for e in r]
    except TypeError as e:
        return r


def all_residues_sans_water(imol):
    return residues_matching_criteria(imol,
                                      lambda chain_id, res_no, ins_code, serial: coot.residue_name(imol, chain_id, res_no, ins_code) != "HOH")

# Return a list of all the residues in the chain
#


def residues_in_chain(imol, chain_id_in):
    """Return a list of all the residues in the chain"""
    return residues_matching_criteria(imol,
                                      lambda chain_id, resno, ins_code, serial: chain_id == chain_id_in)

# Return a list of all residues that have alt confs: where a residue
# is specified thusly: [[chain_id, resno, ins_code], [...] ]
#


def residues_with_alt_confs(imol):

    # return False if there are no atoms with alt-confs, else return
    # a list of the residue's spec [chain_id, res_no, ins_code]
    def alt_func1(chain_id, res_no, ins_code, res_serial_no):
        r = False
        atom_ls = coot.residue_info_py(imol, chain_id, res_no, ins_code)
        for i in range(len(atom_ls)):
            alt_conf_str = atom_ls[i][0][1]
            if alt_conf_str:
                r = True
                break
        return r

    return residues_matching_criteria(imol, lambda chain_id, res_no, ins_code, res_serial_no: alt_func1(chain_id, res_no, ins_code, res_serial_no))

# Return a list of all the altconfs in the residue.
# Typically this will return [""] or ["A", "B"]
#


def residue_alt_confs(imol, chain_id, res_no, ins_code):

    atom_ls = coot.residue_info_py(imol, chain_id, res_no, ins_code)
    alt_confs = []

    if atom_ls:
        for i in range(len(atom_ls)):
            alt_conf_str = atom_ls[i][0][1]
            if not alt_conf_str in alt_confs:
                alt_confs.append(alt_conf_str)
    return alt_confs


# Return a list of all atoms that have zero occupancy: where an atom
# is specified thusly: [[chain_id, resno, ins_code, name, altconf], [...] ]
#
def atoms_with_zero_occ(imol):

    r = []
    for chain_id in chain_ids(imol):
        n_residues = coot.chain_n_residues(chain_id, imol)
        for serial_number in range(n_residues):

            res_name = coot.resname_from_serial_number(
                imol, chain_id, serial_number)
            res_no = coot.seqnum_from_serial_number(imol, chain_id, serial_number)
            ins_code = coot.insertion_code_from_serial_number(
                imol, chain_id, serial_number)
            res_info = coot.residue_info_py(imol, chain_id, res_no, ins_code)
            for atom_info in res_info:
                occ = atom_info[1][0]
                name = atom_info[0][0]
                alt_conf = atom_info[0][1]
                if (occ < 0.01):
                    text = str(imol) + " " + chain_id + " " + \
                        str(res_no) + " " + name
                    if (alt_conf):
                        text += " "
                        text += alt_conf
                    r.append([text, imol, chain_id, res_no,
                              ins_code, name, alt_conf])

    return r

# not to be confused with residue_atom_to_atom_name
# (which uses the output of coot.residue_info_py)
#
# extraction function
#


def atom_spec_to_chain_id(atom_spec):
    """Extract chain ID from atom spec.
    
    Parameters
    ----------
    atom_spec : list
        Atom spec [imol, chain_id, res_no, ins_code, atom_name, alt_conf, ...]
        
    Returns
    -------
    str
        Chain ID
    """
    # atom_spec example ["A", 7, "", " SG ", ""]
    ret = False
    if not atom_spec:
        return False
    if (len(atom_spec) == 5):
        return atom_spec[0]
    elif (len(atom_spec) == 6):
        return atom_spec[1]
    else:
        return False

# extraction function


def atom_spec_to_res_no(atom_spec):
    # atom_spec example ["A", 7, "", " SG ", ""]
    ret = False
    if not atom_spec:
        return False
    if (len(atom_spec) == 5):
        return atom_spec[1]
    elif (len(atom_spec) == 6):
        return atom_spec[2]
    else:
        return False

# extraction function


def atom_spec_to_ins_code(atom_spec):
    # atom_spec example ["A", 7, "", " SG ", ""]
    ret = False
    if not atom_spec:
        return False
    if (len(atom_spec) == 5):
        return atom_spec[2]
    elif (len(atom_spec) == 6):
        return atom_spec[3]
    else:
        return False

# extraction function


def atom_spec_to_atom_name(atom_spec):
    # atom_spec example ["A", 7, "", " SG ", ""]
    ret = False
    if not atom_spec:
        return False
    if (len(atom_spec) == 5):
        return atom_spec[3]
    elif (len(atom_spec) == 6):
        return atom_spec[4]
    else:
        return False

# extraction function


def atom_spec_to_alt_loc(atom_spec):
    # atom_spec example ["A", 7, "", " SG ", ""]
    ret = False
    if not atom_spec:
        return False
    if (len(atom_spec) == 5):
        return atom_spec[4]
    if (len(atom_spec) == 6):
        return atom_spec[5]
    else:
        return False


# simple extraction function
#
def res_spec_to_chain_id(res_spec):
    """simple extraction function"""

    if not res_spec:
        return False
    if (len(res_spec) == 4):
        return res_spec[1]
    if (len(res_spec) == 3):
        return res_spec[0]
    return False

# simple extraction function


def res_spec_to_res_no(res_spec):
    """simple extraction function"""

    if not res_spec:
        return False
    if (len(res_spec) == 4):
        return res_spec[2]
    if (len(res_spec) == 3):
        return res_spec[1]
    return False

# simple extraction function


def res_spec_to_ins_code(res_spec):
    """simple extraction function"""

    if not res_spec:
        return False
    if (len(res_spec) == 4):
        return res_spec[3]
    if (len(res_spec) == 3):
        return res_spec[2]
    return False


# Return False if no atom can be found given the spec else return a list
# consisting of the atom name and alt-conf specifier.
#
# Choose an atom that is called " CA ".  Failing that choose the
# first atom.
#
def residue_spec_to_atom_for_centre(imol, chain_id, res_no, ins_code):

    # residue-info can return False
    atom_ls = coot.residue_info_py(imol, chain_id, res_no, ins_code)

    centre_atom_name_alt_conf = False

    if atom_ls:
        for i in range(len(atom_ls)):
            alt_conf_str = atom_ls[i][0][1]
            atom = atom_ls[i][0][0]
            if (atom == " CA "):
                centre_atom_name_alt_conf = [atom, alt_conf_str]
                break
        if (not centre_atom_name_alt_conf) and (len(atom_ls) > 0):
            # take the first atom
            centre_atom_name_alt_conf = atom_ls[0][0][0:2]

    return centre_atom_name_alt_conf


def set_go_to_atom(res_spec):
    coot.set_go_to_atom_chain_residue_atom_name(
        res_spec_to_chain_id(res_spec),
        res_spec_to_res_no(res_spec),
        " CA ")


def update_go_to_atom_from_current_atom():

    active_atom = active_residue()
    if active_atom:
        imol = active_atom[0]
        chain_id = active_atom[1]
        resno = active_atom[2]
        ins_code = active_atom[3]
        atom_name = active_atom[4]
        alt_conf = active_atom[5]
        go_to_atom_imol_current = coot.go_to_atom_molecule_number()
        coot.set_go_to_atom_molecule(imol)
        # if imol != goto_atom_imol_current
        coot.update_go_to_atom_window_on_other_molecule_chosen(imol)
        coot.set_go_to_atom_chain_residue_atom_name(chain_id, resno, atom_name)


def flip_active_ligand():
    active_atom = active_residue()
    imol = active_atom[0]
    chain_id = active_atom[1]
    resno = active_atom[2]
    ins_code = active_atom[3]
    atom_name = active_atom[4]
    alt_conf = active_atom[5]
    coot.flip_ligand(imol, chain_id, resno)

# Typically one might want to use this on a water, but it deletes the
# nearest CA currently...  Needs a re-think.  Should active-atom just
# return the nearest atom and not be clever about returning a CA.
#


def delete_atom_by_active_residue():

    active_atom = active_residue()
    if active_atom:
        coot.delete_atom(active_atom)

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


def mutate_by_overlap(imol, chain_id_in, resno, tlc):

    # residue is standard residues or phosphorylated version
    #
    # BL says:: maybe this can/should be a global function
    #
    def is_amino_acid(imol, ch_id, res_no):

        aa_list = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLY", "GLU", "GLN",
                   "PHE", "HIS", "ILE", "LEU", "LYS", "MET", "PRO", "SER",
                   "TYR", "THR", "VAL", "TRP", "SEP", "PTR", "TPO", "MSE"]
        rn = coot.residue_name(imol, ch_id, res_no, "")
        if not isinstance(rn, str):
            return False
        else:
            if rn in aa_list:
                return True
        return False

    def is_nucleotide(imol, ch_id, res_no):
        type_list = ["G", "A", "T", "U", "C", "DA", "DG", "DT", "DC"]
        rn = coot.residue_name(imol, ch_id, res_no, "")
        return rn in type_list

    #
    def overlap_by_main_chain(imol_mov, chain_id_mov, res_no_mov, ins_code_mov,
                              imol_ref, chain_id_ref, res_no_ref, ins_code_ref):

        coot.clear_lsq_matches()
        map(lambda atom_name:
            coot.add_lsq_atom_pair_py([chain_id_ref, res_no_ref, ins_code_ref, atom_name, ""],
                                      [chain_id_mov, res_no_mov, ins_code_mov, atom_name, ""]), [" CA ", " N  ", " C  "])
        coot.apply_lsq_matches_py(imol_ref, imol_mov)

    def is_purine(res_name):
        type_list = ["G", "A", "DA", "DG"]
        return res_name in type_list

    def is_pyrimidine(res_name):
        type_list = ["T", "U", "C", "DT", "DC"]
        return res_name in type_list

    def overlap_by_base(imol_mov, chain_id_mov, res_no_mov, ins_code_mov,
                        imol_ref, chain_id_ref, res_no_ref, ins_code_ref):

        coot.clear_lsq_matches()
        rn_1 = coot.residue_name(imol_mov, chain_id_mov, res_no_mov, ins_code_mov)
        rn_2 = coot.residue_name(imol_ref, chain_id_ref, res_no_ref, ins_code_ref)
        purine_set = [" N9 ", " N7 ", " C5 ", " N1 ", " N3 "]
        pyrimidine_set = [" N1 ", " C5 ", " N3 "]
        purine_to_pyrimidine_set = [" N1 ", " C2 ", " N3 "]
        pyrimidine_to_purine_set = [" N9 ", " C4 ", " N5 "]

        atom_list_1 = []
        atom_list_2 = []

        if is_purine(rn_1):
            if is_purine(rn_2):
                atom_list_1 = purine_set
                atom_list_2 = purine_set

        if is_pyrimidine(rn_1):
            if is_pyrimidine(rn_2):
                atom_list_1 = pyrimidine_set
                atom_list_2 = pyrimidine_set

        if is_purine(rn_1):
            if is_pyrimidine(rn_2):
                atom_list_1 = purine_to_pyrimidine_set
                atom_list_2 = pyrimidine_to_purine_set

        if is_pyrimidine(rn_1):
            if is_purine(rn_2):
                atom_list_1 = pyrimidine_to_purine_set
                atom_list_2 = purine_to_pyrimidine_set

        for atom_name_1, atom_name_2 in zip(atom_list_1, atom_list_2):
            atom_spec_r = [chain_id_ref, res_no_ref, ins_code_ref, atom_name_1, ""]
            atom_spec_m = [chain_id_mov, res_no_mov, ins_code_mov, atom_name_2, ""]
            coot.add_lsq_atom_pair_py(atom_spec_r, atom_spec_m)

        coot.apply_lsq_matches_py(imol_ref, imol_mov)
        
    # get_monomer_and_dictionary, now we check to see if we have a
    # molecule already loaded that matches this residue, if we have,
    # then use it.
    #
    def get_monomer_and_dictionary(tlc):

        have_tlc_molecule = False
        for imol in model_molecule_list():
            nc = coot.n_chains(imol)
            if (nc == 1):
                ch_id = chain_id(imol, 0)
                nr = coot.chain_n_residues(ch_id, imol)
                if (nr == 1):
                    rn = coot.resname_from_serial_number(imol, ch_id, 0)
                    if (rn == tlc):
                        have_tlc_molecule = imol
                        break
        have_dict_for_tlc = monomer_restraints(tlc)
        if ((not have_tlc_molecule) or (not have_dict_for_tlc)):
            return coot.get_monomer(tlc)
        else:
            print("we have dict and model for tlc already")
            return have_tlc_molecule

    #
    def mutate_it():
        imol_ligand = get_monomer_and_dictionary(tlc)
        if not valid_model_molecule_qm(imol_ligand):
            s = " Oops.  Failed to get monomer " + str(tlc)
            coot.add_status_bar_text(s)
        else:
            coot.delete_residue_hydrogens(imol_ligand, "A", 1, "", "")
            coot.delete_atom(imol_ligand, "A", 1, "", " OXT", "")

            # 20230207-PE match torsion before main-chain superposition
            if not is_nucleotide(imol_ligand, "A", 1):
                coot.match_ligand_torsions(imol_ligand, imol, chain_id_in, resno)

            if is_amino_acid(imol_ligand, "A", 1) and is_amino_acid(imol, chain_id_in, resno):
                overlap_by_main_chain(imol_ligand, "A", 1, "", imol, chain_id_in, resno, "")
            else:
                overlap_ligands(imol_ligand, imol, chain_id_in, resno)

            # 20200412-PE-merge-complexity:
            if is_nucleotide(imol, chain_id_in, resno):
                if residue_exists_qm(imol, chain_id_in, resno-1, ""):
                    coot.delete_atom(imol_ligand, "A", 1, "", " OP3", "")

            if (is_nucleotide(imol_ligand, "A", 1) and
                is_nucleotide(imol, chain_id_in, resno)):
                overlap_by_base(imol_ligand, "A", 1, "", imol, chain_id_in, resno, "")
            else:
                overlap_ligands(imol_ligand, imol, chain_id_in, resno)

            coot.delete_residue(imol, chain_id_in, resno, "")
            new_chain_id_info = merge_molecules([imol_ligand], imol)
            merge_status = new_chain_id_info[0]

            # merge_status is sometimes a spec, sometimes a chain-id pair
            # BL says:: the rest of it is, merge_status should be 1
            if merge_status == 1:
                new_res_spec = new_chain_id_info[1]
                new_chain_id = residue_spec_to_chain_id(new_res_spec)
                coot.change_residue_number(imol, new_chain_id,
                                      residue_spec_to_res_no(new_res_spec),
                                      residue_spec_to_ins_code(new_res_spec), resno, "")

                if not (new_chain_id == chain_id_in):
                    coot.change_chain_id(imol, new_chain_id, chain_id_in, 1,
                                    residue_spec_to_res_no(new_res_spec),
                                    residue_spec_to_res_no(new_res_spec))

                replacement_state = coot.refinement_immediate_replacement_state()
                imol_map = coot.imol_refinement_map()

                # 20230207-PE doing refinement with a residue that that type non-polymer
                # is bad news. Let's skip that for now

                # coot.set_refinement_immediate_replacement(1)
                # if imol_map == -1:
                #     coot.regularize_zone(imol, chain_id_in, resno, resno, "")
                # else:
                #     spin_atoms = [" P  ", " O1P", " O2P", " O3P"]
                #     phos_dir = {
                #         'PTR': [" CZ ", " OH "],
                #         'SEP': [" CB ", " OG "],
                #         'TPO': [" CB ", " OG1"]}
                #     if tlc in list(phos_dir.keys()):
                #         dir_atoms = phos_dir[tlc]
                #     else:
                #         dir_atoms = False
                #     #refine_zone(imol, chain_id_in, resno, resno, "")
                #     if dir_atoms:
                #         spin_search(imol_map, imol, chain_id_in, resno, "", dir_atoms, spin_atoms)
                #         #refine_zone(imol, chain_id_in, resno, resno, "")
                # coot.accept_regularizement()
                # coot.set_refinement_immediate_replacement(replacement_state)

                coot.set_mol_displayed(imol_ligand, 0)
                coot.set_mol_active(imol_ligand, 0)

            else:
                # guess merge failed!?!
                print("BL WARNING:: merge failed!?")

    # main line
    #
    # First, if there are multiple maps, force the user to choose one,
    # rather than continuing.
    #
    imol_map = coot.imol_refinement_map()
    if (imol_map == -1):
        map_mols = map_molecule_list()
        if len(map_mols) > 1:
            coot.show_select_map_dialog()
            mutate_it()
        else:
            mutate_it()
    else:
        mutate_it()

# A bit of fun
#


def phosphorylate_active_residue():

    def n_active_models():
        n = 0
        for i in model_molecule_list():
            if (coot.mol_is_active(i)):
                n += 1
                return n

    active_atom = active_residue()
    try:
        imol = active_atom[0]
        chain_id = active_atom[1]
        resno = active_atom[2]
        inscode = active_atom[3]
        res_name = coot.residue_name(imol, chain_id, resno, inscode)

        if res_name == 'TYR':
            mutate_by_overlap(imol, chain_id, resno, "PTR")
        elif res_name == 'SER':
            mutate_by_overlap(imol, chain_id, resno, "SEP")
        elif res_name == 'THR':
            mutate_by_overlap(imol, chain_id, resno, "TPO")
        else:
            s = "Can't Phosphorylate residue of type " + res_name
            coot.info_dialog(s)
    except TypeError as e:
        print(e)
        n_active = n_active_models()
        s = 'WARNING:: unable to get active atom\n'
        s += 'Is your molecule active?                  \n'
        s += 'There are ' + str(n_active) + ' active model molecules.'
        coot.info_dialog(s)

# A function for Overlaying ligands.  The transformation is applied
# to all the atoms of the molecule that contains the moving ligand.
#


def overlay_my_ligands(imol_mov, chain_id_mov, resno_mov,
                       imol_ref, chain_id_ref, resno_ref):

    imol_frag = coot.new_molecule_by_atom_selection(imol_mov,
                                               "//" + chain_id_mov +
                                               "/" + str(resno_mov))
    rtop_i = overlap_ligands(imol_frag, imol_ref, chain_id_ref, resno_ref)
    coot.set_mol_displayed(imol_frag, 0)
    transform_coords_molecule(imol_mov, rtop_i[0])


def label_all_CAs(imol):

    for chain_id in chain_ids(imol):
        if not (is_solvent_chain_qm(imol, chain_id)):
            n_residues = coot.chain_n_residues(chain_id, imol)
            for serial_number in number_list(0, n_residues):
                res_name = coot.resname_from_serial_number(
                    imol, chain_id, serial_number)
                res_no = coot.seqnum_from_serial_number(
                    imol, chain_id, serial_number)
                ins_code = coot.insertion_code_from_serial_number(
                    imol, chain_id, serial_number)
                coot.add_atom_label(imol, chain_id, res_no, " CA ")
    coot.graphics_draw()


def label_all_atoms_in_residue(imol, chain_id, resno, inscode):

    atom_list = coot.residue_info_py(imol, chain_id, resno, inscode)
    try:
        for atom_info in atom_list:
            coot.add_atom_label(imol, chain_id, resno, atom_info[0][0])
    except KeyError as e:
        print(e)
    coot.graphics_draw()


def label_all_active_residue_atoms():

    active_atom = active_residue()
    imol = active_atom[0]
    chain_id = active_atom[1]
    resno = active_atom[2]
    inscode = active_atom[3]

    atom_list = coot.residue_info_py(imol, chain_id, resno, inscode)
    if type(atom_list) is ListType:
        for atom_info in atom_list:
            coot.add_atom_label(imol, chain_id, resno, atom_info[0][0])
        coot.graphics_draw()

# Resets alt confs and occupancies of atoms in residue that have
# orphan alt-loc attributes
#


def sanitise_alt_confs(atom_info, atom_ls):

    # return a matching atom (name match) if it exists.  Else return False
    def name_match_qm(atom_1, atom_ls):
        compound_name_1 = atom_1[0]
        atom_name_1 = compound_name_1[0]
        matchers = []
        for atom in atom_ls:
            compound_name_2 = atom[0]
            atom_name_2 = compound_name_2[0]

            if (atom_name_1 == atom_name_2):
                matchers.append(atom)

        if matchers:
            return matchers
        else:
            return False

    # main body
    imol = atom_info[0]
    chain_id = atom_info[1]
    resno = atom_info[2]
    inscode = atom_info[3]
    atom_attribute_settings = []   # add to this

    for atom in atom_ls:
        compound_name = atom[0]
        atom_name = compound_name[0]
        alt_conf = compound_name[1]
        if (alt_conf != ""):
            matchers = name_match_qm(atom, atom_ls)
            if (len(matchers) == 1):
                atom_attribute_settings.append([imol, chain_id, resno, inscode, atom_name,
                                                alt_conf, "alt-conf", ""])
            else:
                atom_attribute_settings.append([imol, chain_id, resno, inscode, atom_name,
                                                alt_conf, "occ", ((shelx_molecule_qm(imol) and 11.0) or (1.0))])
    if (atom_attribute_settings != []):
        set_atom_attributes(atom_attribute_settings)
        if (residue_info_dialog_displayed_qm()):
            coot.residue_info_dialog(imol, chain_id, resno, inscode)

#


def sanitise_alt_confs_in_residue(imol, chain_id, resno, inscode):

    atom_info = [imol, chain_id, resno, inscode, "dummy", "dummy"]
    atom_ls = coot.residue_info_py(imol, chain_id, resno, inscode)
    sanitise_alt_confs(atom_info, atom_ls)

# Resets alt confs and occupancies of atoms in residue that have
# orphan alt-loc attributes.  Use the active-residue.
#


def sanitise_alt_confs_active_residue():
    active_atom = active_residue()
    if active_atom:
        imol = active_atom[0]
        chain_id = active_atom[2]
        resno = active_atom[2]
        inscode = active_atom[2]

        atom_ls = coot.residue_info_py(imol, chain_id, resno, inscode)

        if atom_ls:
            sanitise_alt_confs(active_atom, atom_ls)


def print_molecule_names():

    list(map(lambda molecule_number: printf("    %s    %s\n" % (molecule_number, coot.molecule_name(molecule_number))),
             molecule_number_list()))

# save the dialog positions to the coot_dialog_positions.py file in ./coot-preferences
#


def save_dialog_positions_to_init_file():

    import os
    # return False on failure to find .coot.py (or .coot-preferences)

    def dump_positions_to_file(positions, preferences=False):
        home = 'HOME'
        port = False
        if (os.name == 'nt'):
            home = 'COOT_HOME'
        if preferences:
            init_dir = os.path.join(os.getenv(home),
                                    ".coot-preferences")
            init_file = os.path.join(init_dir,
                                     "saved_dialog_positions.py")
            if (not os.path.isdir(init_dir)):
                return False
            port = open(init_file, 'w')
            print(
                "BL INFO:: writing dialog positions to .coot-preferences/saved_dialog_positions.py")
        else:
            init_file = os.path.join(os.getenv(home),
                                     ".coot.py")
            if (not os.path.isfile(init_file)):
                return False
            port = open(init_file, 'a')
            print("BL INFO:: writing dialog positions to .coot.py")

        if port:

            port.write(
                "# ----------------------------------------------------------")
            port.write(
                "# the following were written by extraction from a state file")
            port.write(
                "# ----------------------------------------------------------")
            for position in positions:
                port.write(position)
            port.write("# -------------------------")
            port.close()
        else:
            print("BL ERROR:: no valid port to write to!")

    # main line
    coot.save_state()

    # FYI, the graphics window is set using
    #
    # coot.set_graphics_window_size(643, 500)
    # coot.set_graphics_window_position(0, 1)
    # They are not dialogs

    state_file = "0-coot.state.py"
    if (not os.path.isfile(state_file)):
        print(
            "Ooops %s does not exist (either guile enabled or problem writing the file" % state_file)
    else:
        port = open("0-coot.state.py", 'r')
        try:
            lines = port.readlines()
        except:
            lines = []
        positions = []
        for line in lines:
            if ("_dialog_position" in line):
                print(" Adding dialog position: ", line)
                positions.append(line)
            elif ("set_go_to_atom_window_position" in line):
                print(" Adding dialog position: ", line)
                positions.append(line)
        port.close()

        write_init_status = dump_positions_to_file(positions)
        if (write_init_status == False):
            # try to write a file to preferences
            dump_positions_to_file(positions, preferences=True)

# saves a string to a file!
# if the string is already present dont do anything
# optional arg: overwrite
#               False - default
#               True  - overwrite file
#


def save_string_to_file(string, filename, overwrite=False):

    #home = 'HOME'
    # if (os.name == 'nt'):
    #    home = 'COOT_HOME'
    init_file = filename
    if (os.path.isfile(init_file)):
        # init file exists
        if (overwrite):
            port = open(init_file, 'w')
            lines = []
        else:
            port = open(init_file, 'r+')
            lines = port.readlines()
    else:
        port = open(init_file, 'w')
        lines = []
    string_written = False
    for line in lines:
        if string in line:
            # line exists, dont write
            string_written = True
            break
    if (not string_written):
        # append the string
        port.write((string + "\n"))
    port.close()


# removes a line containg all strings of the given list from file
def remove_line_containing_from_file(remove_str_ls, filename):

    #home = 'HOME'
    # if (os.name == 'nt'):
    #    home = 'COOT_HOME'
    #init_file = os.path.join(os.getenv(home), ".coot.py")
    init_file = filename
    if (os.path.isfile(init_file)):
        # init file exists
        port = open(init_file, 'r')
        lines = port.readlines()
        port.close()
    else:
        print("BL INFO:: no %s file, so cannot remove line" % init_file)
        lines = []
    if (lines):
        patt = string.join(remove_str_ls, '|')
        re_patt = re.compile(patt)
        tmp_ls = []
        for line in lines:
            result = re_patt.findall(line)
            if (not len(result) == len(remove_str_ls)):
                tmp_ls.append(line)
        port = open(init_file, 'w')
        lines = port.writelines(tmp_ls)
        port.close()

# multiple maps of varying colour from a given map.
#


def multi_chicken(imol, n_colours=False):

    def rotate_colour_map(col, degrees):
        ret = [degrees/360, col[1], col[2] - degrees/360]
        # ??? not sure about 1st element, think mistake in Paul's scheme script
        return ret

    start = 1.0
    stop = 6.0
    initial_colour = [0.2, 0.2, 0.8]
    colour_range = 360

    if valid_map_molecule_qm(imol):
        if not n_colours:
            n_col = 10
        else:
            n_col = n_colours
        sigma = map_sigma(imol)

        print("range n_col returns: ", list(range(n_col)))

        for icol in range(n_col):
            new_map = coot.copy_molecule(imol)
            frac = icol / float(n_col)   # need a float!!
            contour_level_sigma = start + (stop - start) * frac
            coot.set_last_map_contour_level(sigma * contour_level_sigma)
            rotated_colour = rotate_colour_map(initial_colour, colour_range * frac)
            print("#### initial_colour", initial_colour)
            print("#### rotated_colour", rotated_colour)
            coot.set_last_map_colour(*rotated_colour)
    else:
        print("BL INFO:: %s is not valid map" % imol)


# simple enumeration
#
def BALL_AND_STICK(): return 2

# hilight-colour is specified in degrees (round the colour wheel -
# starting at yellow (e.g. 230 is purple))
#


def hilight_binding_site(imol, centre_residue_spec, hilight_colour, radius):

    if (valid_model_molecule_qm(imol)):

        other_residues = coot.residues_near_residue_py(imol, centre_residue_spec, radius)
        atom_sel_str = residue_spec_to_atom_selection_string(centre_residue_spec)

        imol_new = coot.new_molecule_by_atom_selection(imol, atom_sel_str)
        bb_type = 1
        draw_hydrogens_flag = coot.draw_hydrogens_state(imol)

        coot.set_mol_active(imol_new, 0)
        coot.set_show_environment_distances(1)
        coot.set_molecule_bonds_colour_map_rotation(imol_new, hilight_colour)
        coot.additional_representation_by_attributes(imol_new,
                                                centre_residue_spec[0],
                                                centre_residue_spec[1],
                                                centre_residue_spec[1],
                                                centre_residue_spec[2],
                                                BALL_AND_STICK(),
                                                bb_type, 0.14,
                                                draw_hydrogens_flag)

        list(map(lambda spec: coot.additional_representation_by_attributes(imol,
                                                                      spec[0],
                                                                      spec[1],
                                                                      spec[1],
                                                                      spec[2],
                                                                      BALL_AND_STICK(),
                                                                      bb_type, 0.14,
                                                                      draw_hydrogens_flag),
                 other_residues))


highlight_binding_site = hilight_binding_site  # typo?


# Generate restraints from the residue at the centre of the screen
# using PRODRG. Delete hydrogens from the residue because PRODRG has
# anomalous hydrogens.
#
def prodrg_ify(imol, chain_id, res_no, ins_code):

    new_mol = coot.new_molecule_by_atom_selection(imol,
                                             "//" + chain_id + "/" + str(res_no))

    coot.set_mol_active(new_mol, 0)
    coot.set_mol_displayed(new_mol, 0)

    prodrg_dir = "coot-ccp4"
    res_name = coot.residue_name(imol, chain_id, res_no, ins_code)

    if res_name:
        coot.make_directory_maybe(prodrg_dir)
        prodrg_xyzin = os.path.join(prodrg_dir, "prodrg-in.pdb")
        prodrg_xyzout = os.path.join(prodrg_dir, "prodrg-" + res_name + ".pdb")
        prodrg_cif = os.path.join(prodrg_dir, "prodrg-" + res_name + ".cif")
        prodrg_log = os.path.join(prodrg_dir, "prodrg.log")

        coot.delete_residue_hydrogens(new_mol, chain_id, res_no, ins_code, "")
        coot.delete_residue_hydrogens(imol,    chain_id, res_no, ins_code, "")  # otherwise they fly
        coot.write_pdb_file(new_mol, prodrg_xyzin)
        coot.close_molecule(new_mol)
        prodrg_exe = find_exe("cprodrg", "CBIN", "CCP4_BIN", "PATH")
        if not prodrg_exe:
            coot.info_dialog(
                "Cannot find cprodrg, so no prodrg-ifying of ligand possible")
        else:
            print("BL DEBUG:: now run prodrg with", prodrg_exe,
                  "XYZIN",  prodrg_xyzin,
                  "XYZOUT", prodrg_xyzout,
                  "LIBOUT", prodrg_cif)
            status = popen_command(prodrg_exe,
                                   ["XYZIN",  prodrg_xyzin,
                                    "XYZOUT", prodrg_xyzout,
                                    "LIBOUT", prodrg_cif],
                                   ["MINI PREP", "END"],
                                   prodrg_log,
                                   True)
            if status == 0:
                coot.read_cif_dictionary(prodrg_cif)
                imol_new = coot.handle_read_draw_molecule_with_recentre(
                    prodrg_xyzout, 0)
                rn = coot.residue_name(imol, chain_id, res_no, ins_code)
                with_auto_accept([regularize_zone, imol_new, "", 1, 1, ""])
                overlap_ligands(imol_new, imol, chain_id, res_no)
                coot.match_ligand_torsions(imol_new, imol, chain_id, res_no)  # broken?
                overlap_ligands(imol_new, imol, chain_id, res_no)
                coot.set_residue_name(imol_new, "", 1, "", rn)
                coot.change_chain_id(imol_new, "", chain_id, 1, 1, 1)
                coot.renumber_residue_range(imol_new, chain_id, 1, 1, res_no - 1)
                coot.set_mol_displayed(imol_new, 0)
                coot.set_mol_active(imol_new, 0)
                #set_mol_displayed(imol, 0)
                #set_mol_active   (imol, 0)

                # I don't think that replace-fragment is the right
                # function because that does not copy across the hydrogen
                # atoms - and we want those, probably

                # coot.replace_fragment(imol, imol_new,
                #                 "//" + chain_id + "/" + str(res_no))

                imol_replacing = coot.add_ligand_delete_residue_copy_molecule(imol_new, chain_id, res_no, imol, chain_id, res_no)
                col = coot.get_molecule_bonds_colour_map_rotation(imol)
                new_col = col + 5
                coot.set_molecule_bonds_colour_map_rotation(imol_replacing, new_col)
                coot.set_mol_displayed(imol_replacing, 0)
                coot.set_mol_active(imol_replacing, 0)
                coot.graphics_draw()


# ---------- annotations ---------------------

def add_annotation_here(text):
    global annotations
    rc = rotation_centre()
    ann = [text] + rc

    annotations.append(ann)
    coot.place_text(*(ann + [0]))
    coot.graphics_draw()


def add_annotation_at_click(text):
    global pass_text
    pass_text = text

    def add_here(*args):
        global annotations
        global pass_text
        text = pass_text
        # atom_specs for user_defined_clicks have 7 args!
        # includes model number now too!
        # maybe there should be a atom_spec including model no!?
        atom_spec = atom_specs(*args[0][1:7])
        ann = [text] + atom_spec[3:]
        annotations.append(ann)
        coot.place_text(*(ann + [0]))
        coot.graphics_draw()
    coot.user_defined_click_py(1, add_here)


def save_annotations(file_name):
    global annotations
    # if (os.path.isfile(file_name)):
    # remove existing file
    #    print "BL INFO:: overwrite old annotation file", file_name
    #    os.remove(file_name)

    save_string_to_file(str(annotations), file_name, True)


def load_annotations(file_name):
    if (os.path.isfile(file_name)):
        from types import ListType
        port = open(file_name, 'r')
        ls = port.readline()
        port.close()
        ls = ls.rstrip("\n")
        ls = eval(ls)
        if (type(ls) is ListType):
            global annotations
            annotations = ls
            for ann in annotations:
                coot.place_text(*(ann + [0]))
            coot.graphics_draw()


def remove_annotation_here(rad=1.5):
    args = rotation_centre() + [rad]
    handle = coot.text_index_near_position(*args)
    if handle > -1:
        coot.remove_text(handle)


def remove_annotation_at_click(rad=1.5):
    def remove_here(*args):
        # atom_specs for user_defined_clicks have 7 args!
        # includes model number now too!
        # maybe there should be a atom_spec including model no!?
        atom_spec = atom_specs(*args[0][1:7])
        coords = atom_spec[3:]
        handle = coot.text_index_near_position(*(coords + [rad]))
        if handle > -1:
            coot.remove_text(handle)
    coot.user_defined_click_py(1, remove_here)

# ---------- updating ---------------------


global pending_install_in_place
pending_install_in_place = False
# shall this be global? Currently pass use_curl along? FIXME
#global use_curl
#use_curl = False

# Here we construct the url that contains the latest (pre)
# release info adding in "pre-release" if this binary is a
# pre-release.
# args ends up as something like:
# ["xxx/phone_home.py", "pre-release"
# "binary", "Linux-1386-fedora-10-python-gtk2"
# "command-line", "/home/xx/coot/bin/coot"]
#


def make_latest_version_url():
    build_type = coot_sys_build_type()
    # FIXME this is only for biop files !!!
    # what about versions in York?
    host = "http://www.biop.ox.ac.uk/coot/software/binaries/"
    pre_dir = "pre-releases" if pre_release_qm() else "releases"
    if is_windows():
        host = "http://www.ysbl.york.ac.uk/~lohkamp/software/binaries/"
        pre_dir = "nightlies/pre-release" if pre_release_qm() else "stable"
    url = host + \
        pre_dir + \
        "/" + \
        "type-binary-" + \
        build_type + \
        "-latest.txt"
    return url


# get revision number from string
#
# (used in downloading new version)
#
def get_revision_from_string(stri):
    # e.g. str is "coot-0.6-pre-1-revision-2060" (with a newline at the
    # end too).  We want to return 2060 (a number) from here (or False).
    if stri:   # maybe test for string
        lss = stri.rsplit("\n")
        ls = lss[0].split("-")
        try:
            return int(ls[-1])
        except:
            return False
    return False

# for true pythonic url retrievel
# returns url as string or False if error
#


def get_url_as_string(my_url):

    import urllib.request
    import urllib.parse
    import urllib.error
    try:
        s = urllib.request.urlopen(my_url).read()
    except:
        s = False  # not sure if that's the right way
    return s


# first generate a version string with no trailing newline.
# actually removes last char (\n) and everything before c/W.
#
# e.g. input:  "coot-0.6.2-pre-1-revision-2765\n"
#      output: "coot-0.6.2-pre-1-revision-2765"
#
def coot_split_version_string(stri):
    if not is_windows():
        ls2 = stri[stri.find("c"):-1]  # for 'coot' start
    else:
        ls2 = stri[stri.find("W"):-1]  # for 'WinCoot' start
    return ls2

# convert file to string


def file_to_string(file_name):
    if not os.path.isfile(file_name):
        return False
    else:
        fin = open(file_name)
        ret = fin.read()
        fin.close()
        return ret

# If "default.seq" (a simple text file with the sequence (not PIR or
# FASTA)) exists in the current directory, then try to assign it to
# each chain of each molecule.
#
# In the first case the sequence is assigned to the closest match
# (model sequence to target sequence), subsequently only chains
# without a sequence associated with them are candidates for
# matching.  The protein sequence has to have at least 95% sequence
# identity with the target sequence in "default.seq"
#


def load_default_sequence():

    default_seq = "default.seq"
    if os.path.isfile(default_seq):
        s = file_to_string(default_seq)
        coot.align_to_closest_chain(s, 0.95)


def use_curl_status():
    global use_curl
    return use_curl


def set_use_curl(status):
    global use_curl
    use_curl = status


# Invert the chiral centre of the atom we are centred on.
# If not centred on a chiral atom, then give a dialog.
#
# The restraints for this monomer type are copied and renamed
# (changing comp-id, 3-letter-code and name too).  The monomer is
# regularized.  Chiral Hydrogen (if needed) is enabled now in the
# minimizer.
#
# This should almost all be c++ code so that Bernie doesn't have to
# redo it.  This is temporary then
#
# remove 29/11/14
# not needed any more
# def chiral_centre_inverter():
# just to do something. Wait until this is in c++ code...
##    coot.info_dialog("BL waiting for PE to put in C++ code.\n")
# return False

def residue_is_close_to_screen_centre_qm(imol, chain_id, res_no, ins_code):
    def square(x): return x * x
    rc = residue_centre(imol, chain_id, res_no, ins_code)
    if not isinstance(rc, list):
        return False
    else:
        sc = rotation_centre()
        dist_sum_sq = sum(square(rcx-scx) for rcx, scx in zip(rc, sc))
        return dist_sum_sq < 25.

from xml.etree import ElementTree

def fetch_drug_via_wikipedia(drug_name_in):

    def parse_wiki_drug_xml(tree, key):
        query_ele = tree.find("query")
        rev_iter = query_ele.iter("rev")
        for item in rev_iter:
            # print("item", item)
            item_text = item.text
            for idx,line in enumerate(item_text.split("\n")):
                # print("searching for", key, "in", line)
                # try to find the "| DrugBank = DB0001" without finding anything else
                if key in line:
                    if "|" in line:
                        if " = " in line:
                            print("DEBUG:: ", idx, line)
                            parts = line.split()
                            if len(parts) > 2:
                                if len(parts) < 5:
                                    code = parts[-1]
                                    return code
        return False

    def get_redirected_drug_name(xml_tree):
        query_ele = xml_tree.find("query")
        rev_iter = query_ele.iter("rev")
        for item in rev_iter:
            if "#REDIRECT " in item.text:
                lines = item.text.split("\n")
                for line in lines:
                    if "#REDIRECT " in line:
                        parts = line.split()
                        last_word = parts[-1]
                        print("last_word", last_word);
                        p_1 = last_word.strip("[")
                        print("p_1", p_1);
                        p_2 = p_1.strip("]")
                        print("p_2:", p_2);
                        redirected_drug_name = p_2
                        return redirected_drug_name
        return False

    drug_name = drug_name_in.lower()

    def get_xml_from_wikipedia(drug_name):
        # import requests
        url = "http://en.wikipedia.org/w/api.php?format=xml&action=query&titles=" + drug_name + \
            "&prop=revisions&rvprop=content"
        print("DEBUG:: url", url)
        fn = drug_name + ".xml" # put this in a XDG directory.
        state = coot.coot_get_url(url, fn)
        if state == 0:
            # xml_tree = ElementTree.fromstring(xml.content)
            xml_tree = ElementTree.parse(fn)
            return xml_tree
        else:
            return False

    xml_tree = get_xml_from_wikipedia(drug_name)

    redirected_drug_name = get_redirected_drug_name(xml_tree)
    print("redirected_drug_name", redirected_drug_name)
    if redirected_drug_name:
        xml_tree = get_xml_from_wikipedia(redirected_drug_name.lower())

    db_code = parse_wiki_drug_xml(xml_tree, "DrugBank ")
    print("db_code", db_code)
    if db_code:
        db_mol_uri = "https://www.drugbank.ca/structures/small_molecule_drugs/" + str(db_code) + ".mol"
        file_name = str(db_code) + ".mol"
        print("DEBUG:: file_name", file_name)
        print("DEBUG:: db_mol_uri", db_mol_uri)
        coot.coot_get_url(db_mol_uri, file_name)
        return file_name
    print("DEBUG:: returning False")
    return False

def get_SMILES_for_comp_id_from_pdbe(comp_id):

    if not isinstance(comp_id, str):
        return False
    else:
        s = coot.SMILES_for_comp_id(comp_id)
        if isinstance(s, str):
            return s
        else:
            cif_file_name = get_pdbe_cif_for_comp_id(comp_id)
            if isinstance(cif_file_name, str):
                coot.read_cif_dictionary(cif_file_name)
                s = coot.SMILES_for_comp_id(comp_id)
                return s

# return False or a file_name
#


def get_pdbe_cif_for_comp_id(comp_id):
    """return False or a file_name"""

    download_dir = get_directory("coot-download")
    cif_file_name = os.path.join(download_dir,
                                 "PDBe-" + comp_id + ".cif")
    url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/mmcif/" + \
          comp_id + ".cif"

    if os.path.isfile(cif_file_name):
        # try the filesystem cache
        stat_data = os.stat(cif_file_name)
        l = stat_data.st_size
        if (l > 0):
            return cif_file_name
        else:
            # give a dialog, saying that the file will not be
            # overwritten
            msg = cif_file_name + " exists but is empty." + \
                "\nNot overwriting."
            coot.info_dialog(msg)
            return False
    # use network then
    print("INFO:: getting url:", url)
    state = coot.coot_get_url(url, cif_file_name)
    if (state != 0):
        msg = "Problem downloading\n" + url + "\n to file \n" + \
              cif_file_name + "."
        coot.info_dialog(msg)
        return False
    else:
        # it worked!?!
        return cif_file_name

    # something probably went wrong if we got to here
    return False

# # load the redefining functions
# try:
#     load_from_search_load_path("redefine_functions.py")
#     # import redefine_functions
# except:
#     print "load_from_search_load_path() of redefine_functions.py failed"
#     pass

# Add file with filename to preferences directory
#


def file_to_preferences(filename):
    """Copy the file filename from python pkgdatadir directory to
    prefereneces directory.
    """

    import shutil

    coot_python_dir = os.getenv("COOT_PYTHON_DIR")
    if not coot_python_dir:
        coot_python_dir = os.path.join(sys.prefix, 'lib', 'python2.7', 'site-packages', 'coot')
        if not os.path.isdir(coot_python_dir) and is_windows():
            # maybe its old an in another place!?!
            coot_python_dir = os.path.normpath(os.path.join(sys.prefix,
                                                            'lib', 'site-packages', 'coot'))
        else:
            coot_python_dir = os.path.join(
                sys.prefix, 'lib', 'python2.7', 'site-packages', 'coot')

    if not os.path.isdir(coot_python_dir):
        coot.add_status_bar_text("Missing COOT_PYTHON_DIR")
    else:
        ref_py = os.path.join(coot_python_dir, filename)

        if not os.path.exists(ref_py):
            coot.add_status_bar_text("Missing reference template key bindings.")
        else:
            # happy path
            home = os.getenv("HOME")
            if is_windows():
                home = os.getenv("COOT_HOME")
                if not home:
                    # try HOME
                    home = os.getenv("HOME")
                else:
                    # fallback
                    home = os.getenv("USERPROFILE")
            if isinstance(home, str):
                pref_dir = os.path.join(home, ".coot-preferences")
                if not os.path.isdir(pref_dir):
                    coot.make_directory_maybe(pref_dir)
                if not os.path.isdir(pref_dir):
                    coot.add_status_bar_text(
                        "No preferences dir, no keybindings. Sorry")
                else:
                    pref_file = os.path.join(pref_dir, filename)
                    # don't install it if it is already in place.
                    if os.path.isfile(pref_file):
                        s = "keybinding file " + pref_file + \
                            " already exists. Not overwritten."
                        coot.add_status_bar_text(s)
                    else:
                        # check the directory first
                        if not os.path.isdir(pref_dir):
                            coot.make_directory_maybe(pref_dir)
                        shutil.copyfile(ref_py, pref_file)
                        if os.path.isfile(pref_file):
                            exec(compile(open(pref_file, "rb").read(),
                                         pref_file, 'exec'), globals())



# something like this for intermediate atoms also?
#
# n-neighbs is either 0 or 1
#
def rebuild_residues_using_db_loop(imol, middle_residue_spec, n_neighbs):

    global db_loop_preserve_residue_names

    # utility function
    #
    def remove_any_GLY_CBs(imol_db_loop, saved_residue_names, ch_id, resno_low):

        for residue_idx, res_name in enumerate(saved_residue_names):
            if isinstance(res_name, str):
                if res_name == "GLY":
                    res_no_gly = resno_low + residue_idx
                    coot.delete_atom(imol_db_loop, ch_id, res_no_gly, "", " CB ", "")

    # main line

    resno_mid = residue_spec_to_res_no(middle_residue_spec)
    resno_low = resno_mid - n_neighbs
    resno_high = resno_mid + n_neighbs
    r = list(range(resno_mid - 4 - n_neighbs, resno_mid + n_neighbs + 0)) + \
        list(range(resno_mid + 1 + n_neighbs, resno_mid + 3 + n_neighbs))
    ch_id = residue_spec_to_chain_id(middle_residue_spec)
    residue_specs = map(lambda res_no: [ch_id, res_no, ""], r)
    db_loop_preserve_residue_names = True

    loop_mols = protein_db_loops(imol, residue_specs, coot.imol_refinement_map(),
                                 1, db_loop_preserve_residue_names)
    residue_spec_of_residues_to_be_replaced = [middle_residue_spec] \
        if n_neighbs == 0 else map(lambda rno: [ch_id, rno, ""],
                                   range(resno_low, resno_high + 1))

    saved_residue_names = map(lambda rno: residue_spec_to_residue_name(imol, rno),
                              residue_spec_of_residues_to_be_replaced)

    if len(loop_mols) > 0:
        # loop-mols have the correct residue numbering but the wrong chain-id
        # and residue type
        imol_db_loop = loop_mols[1][0]
        tmp_loop_mols = loop_mols[0]
        chain_id_db_loop = chain_id(imol_db_loop, 0)
        if not ch_id == chain_id_db_loop:
            change_chain_id(imol, imol_db_loop, chain_id_db_loop, ch_id, 0, 0, 0)
        selection = "//" + ch_id + "/" + str(resno_low) + "-" + str(resno_high)

        # if the original residue was a GLY, the imol-db-loop might (probably
        # will) have CB. If that is the case, then we should remove the CBs now
        #
        remove_any_GLY_CBs(imol_db_loop, saved_residue_names, ch_id, resno_low)

        # this moves atoms, adds atoms if needed, doesn't change the residue name
        #
        replace_fragment(imol, imol_db_loop, selection)
        # tidy up
        for i in tmp_loop_mols:
            close_molecule(i)

def go_to_box_middle():

    ls = map_molecule_list()
    if ls:
        imol_map = ls[0]
        c = cell(imol_map)
        rc = map(lambda a: a * 0.5, c[:3])
        set_rotation_centre(*rc)

# add terminal residue is the normal thing we do with an aligned
# sequence, but also we can try ton find the residue type of a
# residue in the middle of the chain that is modelled as an ALA, say.
#
# PE comment - not sure why needed.
find_aligned_residue_type = find_terminal_residue_type


def using_gui():
    # we shall see if coot_main_menubar is defined in guile or python
    ret = False
    if coot.coot_has_guile():
        ret = coot.run_scheme_command("(defined? 'coot-main-menubar)")
    if not ret:
        ret = "coot_python" in globals()  # coot_main_menubar is not a global var
    return ret

############################################################################################
# end of Paul's scripting
############################################################################################
#
# some BL functions
#
############################################################################################
#

# for easier switching on of GL lighting on surfaces:
#


def GL_light_on():
    coot.set_do_GL_lighting(1)
    coot.do_GL_lighting_state()

# and to turn it off
#


def GL_light_off():
    coot.set_do_GL_lighting(0)
    coot.do_GL_lighting_state()


# Helper functions to set B-Factors

# set B-factor to bval for molecule imol
#
def set_b_factor_molecule(imol, bval):

    for chain_id in chain_ids(imol):
        start_res = coot.seqnum_from_serial_number(imol, chain_id, 0)
        end_res = coot.seqnum_from_serial_number(
            imol, chain_id, coot.chain_n_residues(chain_id, imol) - 1)
        coot.set_b_factor_residue_range(imol, chain_id, start_res, end_res, bval)

# reset B-factor for molecule imol to default value
#


def reset_b_factor_molecule(imol):

    for chain_id in chain_ids(imol):
        start_res = coot.seqnum_from_serial_number(imol, chain_id, 0)
        end_res = coot.seqnum_from_serial_number(
            imol, chain_id, coot.chain_n_residues(chain_id, imol) - 1)
        coot.set_b_factor_residue_range(
            imol, chain_id, start_res, end_res, coot.default_new_atoms_b_factor())

# reset B-factor for active residue to default value
#


def reset_b_factor_active_residue():

    active_atom = active_residue()

    if not active_atom:
        print("No active atom")
    else:
        imol = active_atom[0]
        chain_id = active_atom[1]
        res_no = active_atom[2]
        ins_code = active_atom[3]
        atom_name = active_atom[4]
        alt_conf = active_atom[5]

        coot.set_b_factor_residue_range(
            imol, chain_id, res_no, res_no, coot.default_new_atoms_b_factor())


# BL module to find exe files
# we need this for popen as it requires the full path of the exe file
# we use arguments and keyword:
#
# program_name	: name of exe to find
#
# args (i.e. path_names) : path name to search (usually "PATH", then maybe CCP4_BIN, ...,
#                         can be a single path as well)
# kwargs        : for some extra argumentsto
#        add_extensions=[list]   pass extra extensions to be tested
#        only_extension=str      use only this one to test
#        no_disk_search=bool     Dont bother searching the disk
#        screen_info=bool        print info etc in console
#
# then we search everywhere
#
# on OS where "which" is available we use this first, rather than
# searching in PATH etc.
#
# returns full path of exe when successful, False otherwise
#
def find_exe(program_name, *args, **kwargs):

    import os
    import string

    global search_disk
    search_disk = None
    info = True  # Yeah... no.
    info = False

    # we shall check for full path names first
    if (os.path.isfile(program_name)):
        return os.path.abspath(program_name)

    # if Unix we use which and python's command module to locate the
    # executable (indepent if PATH was given); commands only available on
    # unix! May use subprocess at some point...
    # if the program is not found with which we use the usual way...
    if (os.name == 'posix'):
        import subprocess
        program_exe = subprocess.getoutput('which ' + program_name)
        if (os.path.isfile(program_exe)):
            return program_exe

    if (len(args) > 0):
        path_ls = args
    else:
        # no extra PATH given, should at least check in this dir
        # and in PATH
        path_ls = ["PATH", os.getcwd()]

    # setting of OS specific path properties
    extensions = []
    drives_ls = ["/"]
    program_name_noext = program_name
    # some windows magic
    if (os.name == 'nt'):
        drives_ls = get_windows_drives()
        program_name_noext = strip_extension(program_name)
        file_ext = file_name_extension(program_name)
        # if extenion is explicitly given - only use this one
        # otherwise try all possible ones on Windows, i.e PATHEXT
        if (file_ext):
            extensions = [file_ext]
        else:
            tmp_ext = os.environ["PATHEXT"].split(os.pathsep)
            # list of extensions (no dot) only
            extensions = [ext[1:] for ext in tmp_ext]

    if "only_extension" in kwargs:
        if kwargs["only_extension"]:
            extensions = kwargs["only_extension"]
    if "add_extension" in kwargs:
        extensions += kwargs["add_extensions"]
    if "screen_info" in kwargs:
        info = kwargs["screen_info"]

    program_names = [program_name_noext]
    if extensions:
        program_names = [program_name_noext + "." + ext for ext in extensions]
        # usually we want to have the one with extension, if there is one
        program_names += [program_name_noext]

    for file_name in program_names:
        # search the extra Paths
        for search_path in path_ls:

            if (os.path.isdir(search_path)):
                # we have a single file name, not environ var
                program_exe = os.path.join(search_path, file_name)
                if (os.path.isfile(program_exe)):
                    if info:
                        print("BL INFO:: We found ", program_exe)
                    return program_exe
            else:
                try:
                    primary_path = os.environ[search_path]
                    for path in string.split(primary_path, os.pathsep):
                        program_exe = os.path.join(path, file_name)
                        if (os.path.isfile(program_exe)):
                            if info:
                                print("BL INFO:: We found ", program_exe)
                            return program_exe
                except:
                    if info:
                        print("BL WARNING:: %s not defined!" % search_path)

    # BL says: before we search everywhere we might want to ask
    # the user if he actually wishes to do so!
    # lets insert a small pygtk dialog and ask!
    # only if graphics
    if "no_disk_search" in kwargs:
        no_search = kwargs["no_disk_search"]
    else:
        no_search = False
    search_disk = False
    if (use_gui_qm and not no_search):
        try:
            # no coot_gui stuff in coot_utils. 20220923-PE
            # search_disk = coot_gui.search_disk_dialog(program_name, path_ls)
            pass
        except NameError as e:
            pass
    if search_disk:
        # search everywhere
        for drive in drives_ls:
            for root, dir, file in os.walk(drive):
                program_exe = os.path.join(root, program_name)
                if (os.path.isfile(program_exe)):
                    return program_exe
    else:
        if info:
            print("BL INFO:: we don't search the whole disk for", program_name_noext)

    if info:
        print("BL WARNING:: We cannot find %s anywhere! Program %s won't run!" % (
            program_name_noext, program_name_noext))
    return False

# for running online docs


def open_url(url):
    import webbrowser

    try:
        webbrowser.open(url, 1, 1)
    except:
        print("BL WARNING:: Cannot open the URL %s in webbrowser %s!" %
              (url, webbrowser.get()))


# to reload modules
def reload_module(name):
    import os
    path = os.getenv('COOT_PYTHON_DIR')
    file = os.path.join(path, name)
    exec(compile(open(file, "rb").read(), file, 'exec'))

# to make print a function:


def printf(*args):
    for arg in args:
        print(arg, end=' ')  # use the right python dumass

# to print elements of a list:


def printl(ls):
    list(map(printf, ls))

# Where cmd is e.g. "bltwish"
#       args is list, e.g. [loggraph, "refmac.log"]
#
# in python < 2.4 (and if no logfile)
#
# Returns the pid or False if failed.
#
# in python >= 2.4 (and with logfile)
#
# Returns the process and the open log file object
#
# uses os.spawn if python version < 2.4 otherwise subprocess
#


def run_concurrently(cmd, args=[], data_list=None, logfile=None, screen_flag=False):
    import sys
    import string
    import os
    import subprocess

    major, minor, micro, releaselevel, serial = sys.version_info

    cmd_execfile = ""
    if not command_in_path_qm(cmd):
        print("WARNING:: run_concurrently(): command ", cmd, " not found in $PATH!")
    else:
        cmd_execfile = find_exe(cmd, "CBIN", "CCP4_BIN", "PATH")

    if cmd_execfile:
        if (major >= 2 and minor >= 4):
            # subprocess
            cmd_args = [cmd_execfile] + args
            log = logfile
            if logfile:
                log = open(logfile, 'w')
            try:
                print("DEBUG:: run_concurrently() log", log)
                process = subprocess.Popen(cmd_args, stdin=subprocess.PIPE, stdout=log)
                print("DEBUG:: run_concurrently() process", process)
                if data_list:
                    for data in data_list:
                        d = data + "\n"
                        print("DEBUG:: run_concurrently(): sending data: ", d.encode())
                        process.stdin.write(d.encode())

                pid = process.pid
                
                print("DEBUG:: run_concurrently() pid", pid)

                if log:
                    return (process, log)
                else:
                    return pid
            except KeyError as e:
                print(e)
                print("WARNING:: run_concurrently(): could not run process with args", cmd_args)
                return False

        else:
            # spawn (old)
            try:
                pid = os.spawnv(os.P_NOWAIT, cmd_execfile, [cmd_execfile] + args)
                return pid
            except:
                print("WARNING:: run_concurrently(): could not run program %s with args %s" % (cmd_execfile, args))
                return False
    else:
        print("WARNING:: run_concurrently(): could not find %s, so not running this program" % cmd)
        return False

# python command to see if we have pygtk available
# return True if availabel otherwise False
#


def coot_has_pygtk():
    import sys
    if ('pygtk' in list(sys.modules.keys())):
        return True
    else:
        return False

# python command to see if we have pygobject available
# return True if availabel otherwise False
#


def coot_has_gobject():
    import sys
    if ('gobject' in list(sys.modules.keys())):
        return True
    else:
        return False


# function to kill a process, given the process pid
# return True if managed to kill the process otherwise false
#
def kill_process(pid):
    import os
    import sys
    import signal

    if (os.name == 'nt'):
        # Windows killing tasks
        try:
            # for now use subprocess.call, maybe can use Popen.kill?!
            major, minor, micro, releaselevel, serial = sys.version_info
            if (major >= 2 and minor >= 4):
                # new style
                import subprocess
                ret = subprocess.call("taskkill /F /PID %i" % pid, shell=True)
            else:
                ret = os.system("taskkill /F /PID %i" % pid)
            if (ret == 0):
                # success
                return True
            else:
                return False
        except:
            return False
    else:
        try:
            os.kill(pid, signal.SIGKILL)
            return True
        except:
            return False


# some example function for the toolbutton
# maybe should go in coot_gui.coot_gui!?
def stereo_mono_toggle():
    display_state = coot.stereo_mode_state()
    if (display_state == 0):
        coot.hardware_stereo_mode()
    else:
        coot.mono_mode()


def side_by_side_stereo_mono_toggle():
    display_state = coot.stereo_mode_state()
    if (display_state == 0):
        coot.side_by_side_stereo_mode(0)
    else:
        coot.mono_mode()


def zalman_stereo_mono_toggle():
    display_state = coot.stereo_mode_state()
    if (display_state == 0):
        coot.zalman_stereo_mode()
    else:
        coot.mono_mode()


def switch_stereo_sides():
    factor = -1. * coot.hardware_stereo_angle_factor_state()
    coot.set_hardware_stereo_angle_factor(factor)


def toggle_full_screen(widget=None):
    """ Toggle between full screen and window mode

    Keyword arguments:
    widget -- can be passed from the toolbutton
    """

    if widget:
        if widget.get_active():
            # the button is toggled on
            coot.full_screen(1)
        else:
            coot.full_screen(0)
    else:
        # no alternative for now (could just go by state and change back and forth)
        print("WARNING:: no widget")


def split_active_water():
    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf]:
        coot.split_water(aa_imol, aa_chain_id, aa_res_no, aa_ins_code)

# helper function to test for a number
# returns True if number, otherwise False
#


def isNumber(num):
    """
    helper function to test for a number
    returns True if number, otherwise False
    """
    import numbers
    if isinstance(num, bool):
        # bool are numbers and give "False" results
        return False
    else:
        return isinstance(num, numbers.Number)


# function to merge multiple solvent chains
#
def merge_solvent_chains(imol):

    # first renumber all water chains
    coot.renumber_waters(imol)
    # collect all solvent chains
    solvent_chains = []
    for chain_id in chain_ids(imol):
        if (is_solvent_chain_qm(imol, chain_id)):
            solvent_chains.append(chain_id)

    # now check for overlapping waters and remove
    # maybe this should rather be done in general when merging molecules
    # as well.
    for chain_id in solvent_chains:
        residue_ls = python_representation(imol, [chain_id])[0][0][1]
        for res in residue_ls:
            res_spec = [chain_id, res[0], res[1]]
            if residue_exists_qm(imol, *res_spec):
                near_residues = coot.residues_near_residue_py(imol, res_spec, 0.05)
                if near_residues:
                    # delete
                    for del_res in near_residues:
                        delete_residue_by_spec(imol, del_res)

    # renumber chains after removal:
    coot.renumber_waters(imol)

    # now merge and renumber
    if (len(solvent_chains) > 1):
        master_chain = solvent_chains[0]
        last_prev_water = coot.chain_n_residues(master_chain, imol)
        for chain_id in solvent_chains[1:]:
            n_residues = coot.chain_n_residues(chain_id, imol)
            coot.renumber_residue_range(imol, chain_id, 1,
                                   n_residues, last_prev_water)
            new_start = last_prev_water + 1
            new_end = last_prev_water + n_residues
            coot.change_chain_id(imol, chain_id, master_chain, 1,
                            new_start, new_end)
            last_prev_water = new_end


# delete this function ASAP.
#
# helper to comvert functions to strings
def cmd2str(*args):
    from types import StringType
    func = args[0]
    if callable(func):
        ret = args[0].__name__
    else:
        ret = str(func)
    ret += "("
    for arg in args[1:]:
        if type(arg) is StringType:
            ret += "\"" + arg + "\""
        else:
            ret += str(arg)
        if not arg == args[-1]:
            ret += ", "
    ret += ")"
    return ret


# Function to change the occupancies of alternative conformations
# alt_conf_list is a list of alt_conf string and occ,
# e.g. [["A", 0.8], ["B", 0.2]]
#
def set_alt_conf_occ(imol, chain_id, res_no, ins_code, alt_conf_list):
    """
    Function to change the occupancies of alternative conformations
    alt_conf_list is a list of alt_conf string and occ,
    e.g. [["A", 0.8], ["B", 0.2]]
    """
    # first check if we have alt confs:
    alt_confs = residue_alt_confs(imol, chain_id, res_no, ins_code)
    if len(alt_confs) > 1:
        atom_ls = coot.residue_info_py(imol, chain_id, res_no, ins_code)
        change_list = []
        for i in range(len(atom_ls)):
            alt_conf_str = atom_ls[i][0][1]
            atom_name = atom_ls[i][0][0]
            for (alt_conf_name, alt_conf_occ) in alt_conf_list:
                if (alt_conf_str == alt_conf_name):
                    change_list.append([imol, chain_id, res_no, ins_code, atom_name,
                                        alt_conf_str, "occ", alt_conf_occ])
        set_atom_attributes(change_list)
    else:
        print("BL INFO:: no alt confs in residue",
              imol, chain_id, res_no, ins_code)


# simplyfy check for windows
def is_windows():
    return os.name == 'nt'

# find all windows drives
#


def get_windows_drives():
    import string
    try:
        from ctypes import windll
        drives = []
        bitmask = windll.kernel32.GetLogicalDrives()
        upper = string.ascii_uppercase
        for letter in upper:
            if bitmask & 1:
                drives.append(letter + ":\\")
            bitmask >>= 1
    except:
        import string
        # poor man's version but simple
        # print "BL INFO:: couldnt import ctypes, using simple version to get drives."
        drives = []
        lower = string.ascii_lowercase
        for letter in lower:
            if os.path.isdir(letter + ":"):
                drives.append(letter + ":\\")

    return drives


# clean up pdb file (imol)
# a wrapper for fix_nomenclature errors, sort chains, residues, merge
# solvent chains, renumber waters, etc.
#
def clean_pdb(imol):
    """clean up pdb file (imol)
    a wrapper for fix_nomenclature errors, sort chains, residues, merge
    solvent chains, renumber waters, etc."""

    coot.fix_nomenclature_errors(imol)
    merge_solvent_chains(imol)
    coot.renumber_waters(imol)
    coot.sort_chains(imol)
    coot.sort_residues(imol)


# acronym
merge_water_chains = merge_solvent_chains

# the python run_thread function if no graphics
# if not use_gui_qm:
#     import threading
#     # this has locked, so that no one else can use it
#     global python_thread_return
#     python_thread_return = False

# function to run a python thread with function using
# args which is a tuple
# optionally pass sleep time in ms (default is 20) - usefull
# for computationally expensive threads which may have run longer
# N.B. requires gobject hence in coot_gui.coot_gui.py
#
# def coot_gui.run_python_thread(function, args):

#     class MyThread(threading.Thread):
#         def __init__(self):
#             threading.Thread.__init__(self)
#         def run(self):
#             global python_thread_return
#             python_return_lock = threading.Lock()
#             python_return_lock.acquire()
#             try:
#                 python_thread_return = function(*args)
#             finally:
#                 python_return_lock.release()

#                 MyThread().start()

enhanced_ligand_coot_qm = coot.enhanced_ligand_coot_p

# Function to hide hydrogens in all molecules
#


def hide_all_hydrogens():
    """This will hide all hydrogens. They are not deleted."""

    for imol in model_molecule_number_list():
        coot.set_draw_hydrogens(imol, 0)

# Function to show all available hydrogens. No generation.
#


def show_all_hydrogens():
    """This will show hydrogens on all models, if available. It wont generate any."""

    for imol in model_molecule_number_list():
        coot.set_draw_hydrogens(imol, 1)

# Duplication of a given residue range (in alt conf of course)
#


def duplicate_residue_range(imol, chain_id, res_no_start, res_no_end,
                            occ_split=0.5):
    """
    This function duplicates the given residue range and makes two
    alternative conformations of it. The occupancies are split 50:50 by
    default.

    Args:
        imol: molecule number
        chain_id: chain
        res_no_start: start of residue range
        res_no_end: end of residue range
    Keyword Args:
        occ_split: alt conformation occupancy for alt conf A (default 0.5)
    """

    # backup current state and turn off backup
    occ_backup = coot.get_add_alt_conf_new_atoms_occupancy()
    split_type_backup = coot.alt_conf_split_type_number()
    inter_state = coot.show_alt_conf_intermediate_atoms_state()
    coot.make_backup(imol)  # do a backup first
    backup_mode = coot.backup_state(imol)
    coot.turn_off_backup(imol)

    # set new state
    occ_a = (1 - occ_split) if (occ_split < 1 and occ_split > 0) else 0.5
    coot.set_add_alt_conf_new_atoms_occupancy(1 - occ_split)
    coot.set_add_alt_conf_split_type_number(1)  # better 2 for range?!
    coot.set_show_alt_conf_intermediate_atoms(1)

    ins_code = ""
    alt_conf = ""  # maybe this and above should be arg!?
    rot_no = 0  # ignored currently anyway
    for res_no in range(res_no_start, res_no_end + 1):
        add_alt_conf(imol, chain_id, res_no, ins_code, alt_conf, rot_no)
        coot.accept_regularizement()

    # set things back
    coot.set_add_alt_conf_new_atoms_occupancy(occ_backup)
    coot.set_add_alt_conf_split_type_number(split_type_backup)
    coot.set_show_alt_conf_intermediate_atoms(inter_state)
    if (backup_mode == 1):
        coot.turn_on_backup(imol)

# Necessary for jligand to find libcheck. Mmmh. Was this required before?!
# does similar things to the ccp4 console batch. Win only
# Shouldnt be necessary since we try to setup CCP4 in the startup script.
# And, $CCP4 probably wont be set any more...
#
def setup_ccp4():
    """This will append ccp4 (e.g. CBIN) to PATH, so that we can find
    CCP4 programs in PATH and not only in CBIN. But only if not already
    in PATH."""

    import os

    # only for win
    if (os.name == 'nt'):
        # first check for CCP4
        ccp4_dir = os.getenv("CCP4")
        if not ccp4_dir:
            print("BL INFO:: no $CCP4 found, wont setup CCP4 for Coot.")
            # done here
        else:
            CCP4 = ccp4_dir
            CCP4_MASTER = os.path.abspath(os.path.join(ccp4_dir, os.pardir))
            # not all required I guess!? They should be set anyway
            ccp4_env_vars = {
                "CCP4_SCR": [r"C:\ccp4temp"],
                "CCP4I_TCLTK": [CCP4_MASTER, "TclTk84", "bin"],
                "CBIN": [CCP4, "bin"],
                "CLIB": [CCP4, "lib"],
                "CLIBD": [CCP4, "lib", "data"],
                "CEXAM": [CCP4, "examples"],
                "CHTML": [CCP4, "html"],
                "CINCL": [CCP4, "include"],
                "CCP4I_TOP": [CCP4, "share", "ccp4i"],
                "CLIBD_MON": [CCP4, "lib", "data", "monomers\\"],
                "MMCIFDIC": [CCP4, "lib", "ccp4", "cif_mmdic.lib"],
                "CRANK": [CCP4, "share", "ccp4i", "crank"],
                "CCP4_OPEN": ["UNKNOWN"],
                "GFORTRAN_UNBUFFERED_PRECONNECTED": ["Y"]
            }
            for env_var in ccp4_env_vars:
                env_dir = os.getenv(env_var)
                if not env_dir:
                    # variable not set or empty, so let do so if exists
                    key = ccp4_env_vars[env_var]
                    if len(key) > 1:
                        if os.path.isdir(env_dir):
                            # have dir so set variable
                            key = ccp4_env_vars[env_var]
                            value = os.path.join(key)
                            # print "BL DEBUG:: set env variable to", env_var, value
                            os.environ[env_var] = value
                    else:
                        value = key[0]
                        # print "BL DEBUG:: set env variable to", env_var, value
                        os.environ[env_var] = value
            # change PATH!?
            # how to do this cleverly?! Insert after first 3 - wincoot path
            # (better to test for PATH where Coot comes from)
            path = os.environ["PATH"]
            path_list = path.split(os.pathsep)
            path_list.insert(3, os.path.join(CCP4, "etc"))
            path_list.insert(3, os.path.join(CCP4, "bin"))
            os.environ["PATH"] = os.pathsep.join(path_list)
            # print "BL DEBUG:: PATH set to", os.environ["PATH"]

# Moved from gui_add_linked_cho.py to make a global function.
#


def delete_residue_by_spec(imol, spec):
    coot.delete_residue(imol,
                   residue_spec_to_chain_id(spec),
                   residue_spec_to_res_no(spec),
                   residue_spec_to_ins_code(spec))


# Required if there is no ccp4 in PATH otherwise, e.g. wont find libcheck
# for jligand
#
# 20180603-PE No. This should not be here. Put it it JLigand setup.
# setup_ccp4()
# we work with globals here as to use the function later and not have to bother
# with globals there any more
# set to True in startup to change
global debug_coot
debug_coot = False


def debug():
    global debug_coot
    return debug_coot

# exchange alternative conformations, i.e. in simplest case swap alt conf A
# with alt conf B. Variable to give list of new alt confs names
#


def rename_alt_confs(imol, chain_id, res_no, ins_code,
                     new_names=["B", "A"],
                     old_names=[]):
    """
    Exchange names for alternative conformations, i.e. in simplest case swap
    alt conf A with alt conf B. Variable to give list of new alt conf names.
    E.g. new_names=["B", "A", "D", "C"], swaps A and B as well as C and D.
    Default is alphabetical order for old_names (not necessarily the order
    which is in the residue).
    """

    import string

    if not old_names:
        old_names = list(string.ascii_uppercase)

    # only works if all new alt confs are single letter
    list_len = [len(x) == 1 for x in new_names]
    if not all(list_len):
        print("BL INFO:: not all new alt confs are single letter, method wont work.")
    else:
        # gogogo

        # make dictionary to map new names to old
        alt_conf_dic = dict((old_alt_conf, new_alt_conf)
                            for (old_alt_conf, new_alt_conf)
                            in zip(old_names, new_names))

        # first check if we have alt confs:
        alt_confs = residue_alt_confs(imol, chain_id, res_no, ins_code)
        if len(alt_confs) > 1:
            atom_ls = coot.residue_info_py(imol, chain_id, res_no, ins_code)
            change_list = []
            for i in range(len(atom_ls)):
                alt_conf_str = atom_ls[i][0][1]
                atom_name = atom_ls[i][0][0]
                if alt_conf_str:
                    # first change to new lower case as not to duplicate alt conf
                    change_list.append([imol, chain_id, res_no, ins_code, atom_name,
                                        alt_conf_str, "alt-conf",
                                        alt_conf_dic[alt_conf_str].lower()])
            set_atom_attributes(change_list)
            # make lower case back to upper (in change_list only)
            for i in range(len(change_list)):
                # use new alt conf!
                change_list[i][5] = change_list[i][7]
                change_list[i][7] = change_list[i][7].upper()
            set_atom_attributes(change_list)
        else:
            print("BL INFO:: no alt confs in residue",
                  imol, chain_id, res_no, ins_code)

# swap A and B
#


def rename_alt_confs_active_residue():
    active_atom = active_residue()
    if active_atom:
        imol = active_atom[0]
        chain_id = active_atom[1]
        resno = active_atom[2]
        inscode = active_atom[3]

        rename_alt_confs(imol, chain_id, resno, inscode)


def write_current_sequence_as_pir(imol, ch_id, file_name):
    coot.print_sequence_chain_general(imol, ch_id, 1, 1, file_name)


def run_clustalw_alignment(imol, ch_id, target_sequence_pir_file):

    current_sequence_pir_file = "current-sequence.pir"
    aligned_sequence_pir_file = "aligned-sequence.pir"
    clustalw2_output_file_name = "clustalw2-output-file.log"

    if os.path.exists("aligned-sequence.pir"):
        os.remove("aligned-sequence.pir")
    if os.path.exists("aligned-sequence.dnd"):
        os.remove("aligned-sequence.dnd")
    if os.path.exists("current-sequence.dnd"):
        os.remove("current-sequence.dnd")

    write_current_sequence_as_pir(imol, ch_id, current_sequence_pir_file)

    data_lines = ["3", "1", target_sequence_pir_file, "2", current_sequence_pir_file,
                  "9", "2", "", "4", "", aligned_sequence_pir_file, "", "x", "", "x"]
    popen_command("clustalw2", [], data_lines, clustalw2_output_file_name, 0)
    coot.associate_pir_alignment_from_file(imol, ch_id, aligned_sequence_pir_file)
    coot.apply_pir_alignment(imol, ch_id)
    coot.simple_fill_partial_residues(imol)


# Util function to pipe Coot C stdout to a file (Note: python stdout doesnt
# touch C stdout, therefore this is needed. Of course could just tee out
# all output, but that may not always be required).
#
# based on https://stackoverflow.com/questions/4675728/redirect-stdout-to-a-file-in-python
#
# use either e.g. with open file fn (fn= open(...)):
#
# with stdout_redirected(fn):
#   do some coot things
#
# or:
# with open('output.txt', 'w') as f, stdout_redirected(f):
#  some coot commands
#


import sys
from contextlib import contextmanager

def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd

# @contextmanager - what's this directive?
def stdout_redirected(to=os.devnull, stdout=None):
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = fileno(stdout)
    # copy stdout_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied: 
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(to), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            #NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied



n_alphafold_cols = 60

def make_alphafold_colours():
    cols = []
    for i in range(n_alphafold_cols):
        f = float(i)/float(n_alphafold_cols)
        r = 0.2 + 0.8 * (1.0 - f)
        g = 1.0 * (0.7 - 1.2 * abs(f - 0.5))
        b = 0.1 + 0.9 * f
        ss = r * r + g * g + 0.9 * b * b
        s = math.sqrt(ss)
        r /= s
        g /= s
        b /= s
        cols.append((i, [r,g,b]))
    return cols

def alphafold_pae_to_colour_index(pae):
    f_col = pae * float(n_alphafold_cols) / 100.0
    col_index = int(f_col)
    return col_index

def make_specs_list(file_name):
    specs_list = []
    if os.path.exists(file_name):
        with open(file_name) as f:
            lines = f.readlines()
            for line in lines:
                if line[:4] == "ATOM":
                    if line[12:16] == " CA ":
                        pae_string = line[60:66]
                        pae = float(pae_string)
                        chain_id = line[21:22]
                        res_no_string = line[22:26]
                        res_no = int(res_no_string)
                        ins_code = "" # :-)
                        spec = [chain_id, res_no, ins_code]
                        ci = alphafold_pae_to_colour_index(pae)
                        specs_list.append((spec, ci))
    return specs_list

def make_selections(file_name):
    specs_list = []
    if os.path.exists(file_name):
        with open(file_name) as f:
            lines = f.readlines()
            for line in lines:
                if line[:4] == "ATOM":
                    if line[12:16] == " CA ":
                        pae_string = line[60:66]
                        pae = float(pae_string)
                        chain_id = line[21:22]
                        res_no_string = line[22:26]
                        res_no = int(res_no_string)
                        ins_code = "" # :-)
                        # spec = [chain_id, res_no, ins_code]
                        cid = "//" + chain_id + "/" + str(res_no)
                        ci = alphafold_pae_to_colour_index(pae)
                        specs_list.append((cid, ci))
    return specs_list

def alphafold_pLDDT_colours(imol):
    file_name = coot.molecule_name(imol)
    colours = make_alphafold_colours()
    coot.set_user_defined_colours_py(colours)
    # specs_list = make_specs_list(file_name)
    selections = make_selections(file_name)
    # coot.set_user_defined_atom_colour_by_residue_py(imol, specs_list)
    coot.set_user_defined_atom_colour_by_selection_py(imol, selections)
    coot.graphics_to_user_defined_atom_colours_representation(imol)


def write_current_sequence_as_pir(imol, ch_id, file_name):
    coot.print_sequence_chain_general(imol, ch_id, 1, 1, file_name)

def run_clustalw_alignment(imol, ch_id, target_sequence_pir_file):

    def get_clustalw2_command():
        clustalw2_command = "clustalw2"
        ccp4_libexec_path = os.path.join(os.getenv("CCP4"), "libexec") \
            if os.getenv("CCP4") else ""
        return find_exe(clustalw2_command, ["PATH", ccp4_libexec_path, "CBIN"])

    current_sequence_pir_file = "current-sequence.pir"
    aligned_sequence_pir_file = "aligned-sequence.pir"
    clustalw2_output_file_name = "clustalw2-output-file.log"

    clustalw2_command = get_clustalw2_command()
    if not clustalw2_command:
        print("ERROR:: No clustalw2 command")
        return False
    else:

        if os.path.exists("aligned-sequence.pir"):
            os.remove("aligned-sequence.pir")
        if os.path.exists("aligned-sequence.dnd"):
            os.remove("aligned-sequence.dnd")
        if os.path.exists("current-sequence.dnd"):
            os.remove("current-sequence.dnd")

        write_current_sequence_as_pir(imol, ch_id, current_sequence_pir_file)

        data_lines = ["3", "1", target_sequence_pir_file, "2", current_sequence_pir_file,
                      "9", "2", "", "4", "", aligned_sequence_pir_file, "", "x", "", "x"]
        popen_command("clustalw2", [], data_lines, clustalw2_output_file_name, 0)
        associate_pir_alignment_from_file(imol, ch_id, aligned_sequence_pir_file)
        apply_pir_alignment(imol, ch_id)
        simple_fill_partial_residues(imol)
        resolve_clashing_sidechains_by_deletion(imol)

# Backrub rotamers for chain. After alignment mutation we should run this.
#
# this function cannot be in fitting.py because fitting.py depends on coot_gui
#
def backrub_rotamers_for_chain(imol, ch_id):

    """Backrub rotamers for chain. After alignment mutation we should run this."""

    coot.set_rotamer_search_mode(coot.ROTAMERSEARCHLOWRES)
    coot.make_backup(imol)

    with NoBackups(imol):
        n_times = 2
        imol_map = coot.imol_refinement_map()
        if valid_map_molecule_qm(imol_map):
            n_res = coot.chain_n_residues(ch_id, imol)
            for i_round in range(n_times):
                for serial_number in range(n_res):
                    res_name = coot.resname_from_serial_number(imol, ch_id, serial_number)
                    res_no = coot.seqnum_from_serial_number(imol, ch_id, serial_number)
                    ins_code = coot.insertion_code_from_serial_number(imol, ch_id, serial_number)
                    if isinstance(ins_code, str):   # valid residue check :-)
                        if not res_name == "HOH":
                            coot.auto_fit_best_rotamer(imol, ch_id, res_no, ins_code, "", imol_map, 1, 0.1)




####### Back to Paul's scripting.
####### This needs to follow find_exe

# if you don't have mogul, set this to False
global use_mogul
use_mogul = True
if (not command_in_path_qm("mogul")):
    use_mogul = False

