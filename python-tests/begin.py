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

print "==============================================================="
print "==================== Testing =================================="
print "==============================================================="

import unittest, os
import inspect

global have_test_skip
global skipped_tests

have_test_skip = False
skipped_tests  = []

if ('skip' in dir(unittest.TestCase)):
    have_test_skip = True
else:
    print "WARNING:: unittest skip not avaliable!!!!!!"

home = os.getenv('HOME')
if ((not home) and (os.name == 'nt')):
    home = os.getenv('COOT_HOME')
if (not home):
    # badness we dont have a home dir
    print "ERROR:: Cannot find a HOME directory"


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
    ret = read_pdb(os.path.join(unittest_data_dir, file_name))
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
    ret = math.sqrt(sum(map(square, map(sub, pos_1, pos_2))))
    return ret

pos_diff = bond_length

def bond_length_from_atoms(atom_1, atom_2):
    from types import ListType
    if (type(atom_1) is not ListType):
        print "   WARNING:: bond_length_from_atoms: atom_1 not a list:", atom_1
        return False
    elif (type(atom_2) is not ListType):
        print "   WARNING:: bond_length_from_atoms: atom_2 not a list:", atom_2
        return False
    else:
        return bond_length(atom_1[2],
                           atom_2[2])

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
            test.failUnlessAlmostEqual(occ, 11.0, 1, "  bad occupancy in SHELXL molecule %s" %atom)

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
    from types import ListType, StringType, FloatType, IntType, BooleanType
    if ((not atom_1) or (not atom_2)):
        return False     # no matching residues
    if (len(atom_1) != len(atom_2)):
        return False     # comparing differnt list sizes -> not equal anyway
    for i in range(len(atom_1)):
        if (type(atom_1[i]) is ListType):
            if (atoms_match_qm(atom_1[i], atom_2[i])):
                pass
            else:
                return False
        elif (type(atom_1[i]) is StringType):
            if (atom_1[i] == atom_2[i]):
                pass
            else:
                return False
        elif (type(atom_1[i]) is FloatType):
            if (abs(atom_1[i] - atom_2[i]) < 0.01):
                pass
            else:
                return False
        elif (type(atom_1[i]) is IntType):
            if (atom_1[i] == atom_2[i]):
                pass
            else:
                return False
        elif (type(atom_1[i]) is BooleanType):
            if (atom_1[i] == atom_2[i]):
                pass
            else:
                return False
    return True


# transform
# Just a (eg 3x3) (no vector) matrix
# The matrix does not have to be symmetric.
#
def transpose_mat(mat, defval=None):
    if not mat:
        return []
    return map(lambda *row: [elem or defval for elem in row], *mat)

# What is the distance atom-1 to atom-2? 
# return False on not able to calculate
def atom_distance(atom_1, atom_2):
    def square(x):
        return x*x
    def sub(x, y):
        return x-y
    import math
    ret = math.sqrt(sum(map(square, map(sub, atom_1[2], atom_2[2]))))
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
        return residue_name(imol, *spec[1:4])
    return False

# return residue specs of all residues that are type residue-type
#
def get_residues_in_molecule_of_type(imol, residue_type):

    return filter(lambda x: residue_name_from_spec(imol, x) == residue_type,
                  fit_protein_make_specs(imol, 'all-chains'))

# This takes 4 member specs, return True or False
#
def spec_match_qm(spec_1, spec_2):

    from types import StringType
    from types import IntType
    # first test if specs are of length 4.
    #
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

    if not (type(chain_id_1) is StringType):
        return False
    if not (chain_id_1 == chain_id_2):
        return False
    if not (type(res_no_1) is IntType):
        return False
    if not (res_no_1 == res_no_2):
        return False
    if not (type(ins_code_1) is StringType):
        return False
    if not (ins_code_1 == ins_code_2):
        return False
    return True
    
set_console_display_commands_hilights(0, 0, 0)


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
        
unittest.TestCase.skip_test = skip_test

###############################
# NOW DEFINE THE ACTUAL TESTS
###############################


def get_this_dir():
    pass

test_file_list = ["01_pdb_mtz.py",
                  "02_shelx.py",
                  "03_ligand.py",
                  "04_cootaneer.py",
                  "05_rna_ghosts.py",
                  "06_ssm.py",
                  "07_ncs.py",
                  "08_utils.py",
                  "09_internal.py"]


# get directory of this file and execute tests found in this dir
fn = inspect.getfile(get_this_dir)
current_dir = os.path.dirname(fn)

for test_file in test_file_list:
    load_file = os.path.join(current_dir, test_file)
    load_file = os.path.normpath(load_file)
    if (os.path.isfile(load_file)):
        execfile(load_file, globals())

test_list = [PdbMtzTestFunctions, ShelxTestFunctions,
             LigandTestFunctions, CootaneerTestFunctions,
             RnaGhostsTestFunctions, SsmTestFunctions,
             NcsTestFunctions, UtilTestFunctions,
             InternalTestFunctions]

suite = unittest.TestSuite()
for test in test_list:
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(test))


# returns a list of all tests, in pairs containing the test class and description
# i.e. [[__main__.PdbMtzTestFunctions.test00_0, "Post Go To Atom no molecule"],
#       [__main__.PdbMtzTestFunctions.test01_0, "short discr Close bad molecule"]]
def list_of_all_tests():
    ret =[]
    for test in suite:
        ret.append([test, test.shortDescription(), test.id()])
    return ret
    
global unittest_output
unittest_output = False

# class to write output of unittest into a 'memory file' (unittest_output)
# as well as to sys.stdout
class StreamIO:
        
    def __init__(self, etxra, src=sys.stderr, dst=sys.stdout):
        import StringIO
        global unittest_output
        unittest_output = StringIO.StringIO()
        self.src = src
        self.dst = dst
        self.extra = unittest_output

    def write(self, msg):
        #self.src.write(msg)
        self.extra.write(msg)
        self.dst.write(msg)

    def flush(self):
        pass

class StreamIOnew:

    def __init__(self, etxra, src=sys.stderr, dst=sys.stdout):
        import io
        import StringIO
        global unittest_output
        unittest_output = StringIO.StringIO()
#        unittest_output = io.StringIO()
        self.src = src
        self.dst = dst
        self.extra = unittest_output

    def write(self, msg):
        #msg = msg.decode()
        self.extra.write(msg)
        self.dst.write(msg)

    def flush(self):
        pass

# function to run one test
def run_one_test(no_of_test):
    result = unittest.TextTestRunner(verbosity=2).run(list_of_all_tests()[no_of_test][0])
    return result

# function to run one test set
def run_test_set(no_of_test_set):
    print "BL DEBUG:: should run ", test_list[no_of_test_set]
    if (no_of_test_set > 0):
        read_imol_rnase_and_map()
    set = unittest.TestSuite()
    set.addTests(unittest.TestLoader().loadTestsFromTestCase(test_list[no_of_test_set]))
    result = unittest.TextTestRunner(verbosity=2).run(set)
    return result

# read imol_rnase and map if not there?!
def read_imol_rnase_and_map():
    global imol_rnase
    global imol_rnase_map

    if (imol_rnase < 0):
        # not set up-> read, i.e. run test02 from 01
        imol = read_pdb(rnase_pdb())
        imol_rnase = imol

    if (imol_rnase_map < 0):
        imol_map = make_and_draw_map(rnase_mtz(), "FWT","PHWT","",0,0)
        set_imol_refinement_map(imol_map)
        imol_rnase_map = imol_map
  
