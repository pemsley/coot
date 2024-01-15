
import sys
import unittest
import coot
import coot_utils
import coot_testing_utils

from TestPdbMtzFunctions    import *
from TestShelxFunctions     import *
from TestLigandFunctions    import *
from TestCootaneerFunctions import *
from TestRNAGhostsFunctions import *
from TestSSMFunctions       import *
from TestNCSFunctions       import *
from TestUtilsFunctions     import *
from TestInternalFunctions  import *

# class to write output of unittest into a 'memory file' (unittest_output)
# as well as to sys.stdout
class StreamIO:
        
    def __init__(self, etxra, src=sys.stderr, dst=sys.stdout):
        import io
        global unittest_output
        unittest_output = io.StringIO()
        self.src = src
        self.dst = dst
        self.extra = unittest_output

    def write(self, msg):
        #self.src.write(msg)
        self.extra.write(msg)
        self.dst.write(msg)

    def flush(self):
        pass


suite = unittest.TestSuite()
test_list = [TestPdbMtzFunctions, TestShelxFunctions, TestLigandFunctions, TestCootaneerFunctions,
             TestRNAGhostsFunctions, TestSSMFunctions, TestNCSFunctions, TestUtilsFunctions, TestInternalFunctions]
# test_list = [TestLigandFunctions]

for test in test_list:
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(test))

log = StreamIO(sys.stderr, sys.stdout)

result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
