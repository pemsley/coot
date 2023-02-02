#    begin.py
#    Copyright (C) 2008  Bernhard Lohkamp, The University of York
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

print("===============================================================")
print("==================== Testing ==================================")
print("===============================================================")

import sys
import unittest, os
import coot
import fitting
import inspect

global have_test_skip
global skipped_tests

have_test_skip = False
skipped_tests  = []

print(dir(unittest.TestCase))

if ('skipTest' in dir(unittest.TestCase)):
    have_test_skip = True
else:
    print("WARNING:: unittest skip not avaliable!!!!!!")

home = os.getenv('HOME')
if ((not home) and (os.name == 'nt')):
    home = os.getenv('COOT_HOME')
if (not home):
    # badness we dont have a home dir
    print("ERROR:: Cannot find a HOME directory")


def get_unittest_data_dir():
    d = os.getenv("COOT_TEST_DATA_DIR")
    if (not d):
        d = os.path.normpath(os.path.join(home, "data", "greg-data"))
    return d

global unittest_data_dir
unittest_data_dir = get_unittest_data_dir()   # too lazy to convert all
                                              # occurences of unittest_data_dir

################
# some functions
################

def unittest_pdb(file_name):
    ret = coot.read_pdb(os.path.join(unittest_data_dir, file_name))
    return ret

# use functions rather than variables which would have to be globals
def rnase_pdb():
    return os.path.join(unittest_data_dir, "tutorial-modern.pdb")

def rnase_mtz():
    return os.path.join(unittest_data_dir, "rnasa-1.8-all_refmac1.mtz")

def rotate_n_frames(n):
    rotate_speed = 1
    return int(rotate_speed * n)

# whatever this shall do
#def chains_in_order_qm(ls):

def bond_length(pos_1, pos_2):
    def square(x):
        return x*x
    def sub(x, y):
        return x-y

    import math
    ret = math.sqrt(sum(map(square, list(map(sub, pos_1, pos_2)))))
    return ret

pos_diff = bond_length

def bond_length_from_atoms(atom_1, atom_2):
    # # from types import ListType ; 20230202-PE old
    # if (type(atom_1) is not ListType):
    #     print("   WARNING:: bond_length_from_atoms: atom_1 not a list:", atom_1)
    #     return False
    # elif (type(atom_2) is not ListType):
    #     print("   WARNING:: bond_length_from_atoms: atom_2 not a list:", atom_2)
    #     return False
    # else:
    #     return bond_length(atom_1[2], atom_2[2])

    try:
        return bond_length(atom_1[2], atom_2[2])

    except KeyError as e:
        print(e)
        return False


def bond_length_within_tolerance_qm(atom_1, atom_2, ideal_length, tolerance):

    # BL says:: this should (as in greg) probably throw an error, but we dont have the
    # test object here, so cannot. Buh. Should pass self?!
    if (not atom_1):
        return False
    if (not atom_2):
        return False
    b = bond_length_from_atoms(atom_1, atom_2)
    return abs(b - ideal_length) < tolerance

def shelx_waters_all_good_occ_qm(test, imol_insulin_res):

    chain_id = water_chain(imol_insulin_res)
    n_residues = chain_n_residues(chain_id, imol_insulin_res)
    serial_number = n_residues - 1
    res_name = resname_from_serial_number(imol_insulin_res, chain_id, serial_number)
    res_no   = seqnum_from_serial_number (imol_insulin_res, chain_id, serial_number)
    ins_code = insertion_code_from_serial_number(imol_insulin_res, chain_id, serial_number)

    if (res_name == "HOH"):
        atom_list = residue_info(imol_insulin_res, chain_id, res_no, ins_code)
        for atom in atom_list:
            occ = atom[1][0]
            test.assertAlmostEqual(occ, 11.0, 1, "  bad occupancy in SHELXL molecule %s" %atom)

#  return restraints without the given bond restraints or
# False if no restraints given
def strip_bond_from_restraints(atom_pair, restraints):

    import copy
    restr_dict = copy.deepcopy(restraints)
    if (restr_dict == []):
        return []
    else:
        if not restr_dict:
            return False
        else:
            # make a copy of restraints, so that we dont overwrite
            # the original dictionary
            bonds = restr_dict["_chem_comp_bond"]
            if bonds:
                for i, bond in enumerate(bonds):
                    atom1 = bond[0]
                    atom2 = bond[1]
                    if (atom1 in atom_pair and atom2 in atom_pair):
                        del restr_dict["_chem_comp_bond"][i]
                        return restr_dict
                return restr_dict

# are the attributes of atom-1 the same as atom-2? (given we test for
# failUnless/IfAlmostEqual)
#
def atoms_match_qm(atom_1, atom_2):
    # from types import ListType, StringType, FloatType, IntType, BooleanType
    # if ((not atom_1) or (not atom_2)):
    #     return False     # no matching residues
    # if (len(atom_1) != len(atom_2)):
    #     return False     # comparing differnt list sizes -> not equal anyway
    # for i in range(len(atom_1)):
    #     if (type(atom_1[i]) is ListType):
    #         if (atoms_match_qm(atom_1[i], atom_2[i])):
    #             pass
    #         else:
    #             return False
    #     elif (type(atom_1[i]) is StringType):
    #         if (atom_1[i] == atom_2[i]):
    #             pass
    #         else:
    #             return False
    #     elif (type(atom_1[i]) is FloatType):
    #         if (abs(atom_1[i] - atom_2[i]) < 0.01):
    #             pass
    #         else:
    #             return False
    #     elif (type(atom_1[i]) is IntType):
    #         if (atom_1[i] == atom_2[i]):
    #             pass
    #         else:
    #             return False
    #     elif (type(atom_1[i]) is BooleanType):
    #         if (atom_1[i] == atom_2[i]):
    #             pass
    #         else:
    #             return False
    # return True

    try:
        if len(atom_1) == len(atom_2):
            if atom_1[0] == atom_2[0]:
                if atom_1[1] == atom_2[1]:
                    if atom_1[2] == atom_2[2]:
                        if atom_1[3] == atom_2[3]:
                            print("debug:: atoms_match_qm: a match:", atom_1, atom_2)
                            return True
        else:
            return False
    except KeyError as e:
        print(e)
        return False


# transform
# Just a (eg 3x3) (no vector) matrix
# The matrix does not have to be symmetric.
#
def transpose_mat(mat, defval=None):
    if not mat:
        return []
    return list(map(lambda *row: [elem or defval for elem in row], *mat))

# What is the distance atom-1 to atom-2? 
# return False on not able to calculate
def atom_distance(atom_1, atom_2):
    def square(x):
        return x*x
    def sub(x, y):
        return x-y
    import math
    ret = math.sqrt(sum(map(square, list(map(sub, atom_1[2], atom_2[2])))))
    return ret

# a function from 04_cootaneering:
#
# note that things could go wrong if there is a mising EOL (not tested)
# not sure about this in python
#
def file2string(rnase_pir):
    fin = open(rnase_pir, 'r')
    s = fin.read()
    fin.close()
    return s

def atoms_have_correct_seg_id_qm(atoms, seg_id):

    for atom in atoms:
        atom_seg_id = atom[1][3]
        if (not atom_seg_id == seg_id):
            return False
    return True

# return the residue name given a 4 element spec
#
def residue_name_from_spec(imol, spec):
    
    if len(spec) == 4:
        return coot.residue_name_py(imol, *spec[1:4])
    return False

# return residue specs of all residues that are type residue-type
#
def get_residues_in_molecule_of_type(imol, residue_type):

    return [x for x in fitting.fit_protein_make_specs(imol, 'all-chains') if residue_name_from_spec(imol, x) == residue_type]

# This takes 4 member specs, return True or False
# 20230202-PE I don't like that idea, these days.
#
def spec_match_qm(spec_1, spec_2):

    # from types import StringType
    # from types import IntType
    # # first test if specs are of length 4.
    # #
    # if not (len(spec_1) == 4):
    #     return False
    # if not (len(spec_2) == 4):
    #     return False
    # chain_id_1 = spec_1[1]
    # chain_id_2 = spec_2[1]
    # res_no_1   = spec_1[2]
    # res_no_2   = spec_2[2]
    # ins_code_1 = spec_1[3]
    # ins_code_2 = spec_2[3]

    # if not (type(chain_id_1) is StringType):
    #     return False
    # if not (chain_id_1 == chain_id_2):
    #     return False
    # if not (type(res_no_1) is IntType):
    #     return False
    # if not (res_no_1 == res_no_2):
    #     return False
    # if not (type(ins_code_1) is StringType):
    #     return False
    # if not (ins_code_1 == ins_code_2):
    #     return False
    # return True

    if not (len(spec_1) == 4):
        return False
    if not (len(spec_2) == 4):
        return False

    chain_id_1 = spec_1[1]
    chain_id_2 = spec_2[1]
    res_no_1   = spec_1[2]
    res_no_2   = spec_2[2]
    ins_code_1 = spec_1[3]
    ins_code_2 = spec_2[3]

    if chain_id_1 == chain_id_2:
        if res_no_1 == res_no_2:
            if ins_code_1 == ins_code_2:
                return True

    return False

coot.set_console_display_commands_hilights(0, 0, 0)


#################################
# SOME MORE FUNCTIONS ONLY PYTHON
#################################

def skip_test(self, skip_query, skip_msg):
    global have_test_skip
    global skipped_tests

    if have_test_skip:
        self.skipIf(skip_query, skip_msg)
    else:
        if skip_query:
            skipped_tests.append(self.shortDescription())
            # no way so far to exit the test from here
            # either work with return value or Exceptions
            # former for now
            return True
        
# unittest.TestCase.skip_test = skip_test

