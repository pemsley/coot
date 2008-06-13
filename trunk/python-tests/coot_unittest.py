print "==============================================================="
print "==================== Testing =================================="
print "==============================================================="

import unittest, os
import inspect

home = os.getenv('HOME')
if ((not home) and (os.name == 'nt')):
    home = os.getenv('COOT_HOME')
else:
    # badness we dont have a home dir
    print "ERROR:: Cannot find a HOME directory"

unittest_data_dir = os.path.normpath(os.path.join(home, "data", "greg-data"))

def unittest_pdb(file_name):
    ret = read_pdb(os.path.join(unittest_data_dir, file_name))
    return ret

def rotate_n_frames(n):
    rotate_speed = 1
    return int(rotate_speed * n)

def bond_length(pos_1, pos_2):
    def square(x):
        return x*x
    def sub(x, y):
        return x-y

    import math
    ret = math.sqrt(sum(map(square, map(sub, pos_1, pos_2))))
    return ret

def bond_length_within_tolerance_qm(atom_1, atom_2, ideal_length, tolerance):

    if (not atom_1):
        return False
    if (not atom_2):
        return False
    pos_1 = atom_1[2]
    pos_2 = atom_2[2]
    b = bond_length(pos_1, pos_2)
    return abs(b - ideal_length) < tolerance

def get_atom(imol, chain_id, resno, atom_name):

    def get_atom_from_res(atom_name, residue_atoms):
        for residue_atom in residue_atoms:
            if (residue_atom[0][0] == atom_name):
                return residue_atom
        print "BL WARNING:: no atom name %s found in residue" %atom_name
        return False # no residue name found
    
    res_info = residue_info(imol, chain_id, resno, "")
    if (not res_info):
        return False
    else:
        ret = get_atom_from_res(atom_name, res_info)
        return ret

#  return restraints without the given bond restraints or False if no restraints given
def strip_bond_from_restraints(atom_pair, restraints):

    import copy
    restr_dict = copy.deepcopy(restraints)
    if restr_dict:
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
        return False
    return False
                    

#test_file_list = ["02_shelx.py"]
test_file_list = ["01_pdb_mtz.py",
                  "02_shelx.py",
                  "03_ligand.py",
                  "04_cootaneer.py",
                  "05_rna_ghosts.py",
                  "06_ssm.py",
                  "07_ncs.py",
                  "08_utils.py"]

# get directory of this file and execute tests found in this dir
fn = inspect.getfile(rotate_n_frames)
current_dir = os.path.dirname(fn)

for test_file in test_file_list:
    load_file = os.path.join(current_dir, test_file)
    if (os.path.isfile(load_file)):
        execfile(load_file)

#test_list = [ShelxTestFunctions]
test_list = [PdbMtzTestFunctions, ShelxTestFunctions,
             LigandTestFunctions, CootaneerTestFunctions,
             RnaGhostsTestFunctions, SsmTestFunctions,
             NcsTestFunctions, UtilTestFunctions]

suite = unittest.TestSuite()
for test in test_list:
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(test))

result = unittest.TextTestRunner(verbosity=2).run(suite)

if (result.wasSuccessful()):
    coot_real_exit(0)
else:
    coot_real_exit(1)

