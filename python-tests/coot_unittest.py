print "==============================================================="
print "==================== Testing =================================="
print "==============================================================="

import unittest, os

unittest_data_dir = os.path.join(os.getenv('HOME'), "data", "greg-data")

def rotate_n_frames(n):
    rotate_speed = 1
    return int(rotate_speed * n)

test_file_list = ["01-pdb+mtz.py", "03-ligand.py"]

for test_file in test_file_list:
    execfile(test_file)

test_list = [PdbMtzTestFunctions, LigandTestFunctions]

suite = unittest.TestSuite()
for test in test_list:
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(test))

unittest.TextTestRunner(verbosity=2).run(suite)
