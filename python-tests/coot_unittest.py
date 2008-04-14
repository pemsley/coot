print "==============================================================="
print "==================== Testing =================================="
print "==============================================================="

import unittest, os
import inspect

unittest_data_dir = os.path.normpath(os.path.join(os.getenv('HOME'), "data", "greg-data"))

def rotate_n_frames(n):
    rotate_speed = 1
    return int(rotate_speed * n)

#test_file_list = ["02_shelx.py"]
test_file_list = ["01_pdb_mtz.py",
                  "02_shelx.py",
                  "03_ligand.py",
                  "04_cootaneer.py",
                  "05_rna_ghosts.py",
                  "06_ssm.py",
                  "07_ncs.py"]

# get directory of this file and execute tests found in this dir
fn = inspect.getfile(rotate_n_frames)
current_dir = os.path.dirname(fn)

for test_file in test_file_list:
    execfile(os.path.join(current_dir,test_file))

#test_list = [ShelxTestFunctions]
test_list = [PdbMtzTestFunctions, ShelxTestFunctions,
             LigandTestFunctions, CootaneerTestFunctions,
             RnaGhostsTestFunctions, SsmTestFunctions,
             NcsTestFunctions]

suite = unittest.TestSuite()
for test in test_list:
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(test))

result = unittest.TextTestRunner(verbosity=2).run(suite)

if (result.wasSuccessful()):
    coot_real_exit(0)
else:
    coot_real_exit(1)

