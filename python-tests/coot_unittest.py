#    coot_unittest.py
#    Copyright (C) 2007  Bernhard Lohkamp
#    Copyright (C) 2007, 2008  Bernhard Lohkamp, The University of York
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

import unittest, os
import inspect
import sys
import coot
import coot_utils

def get_this_dir():
    pass


# get directory of this file and execute tests found in this dir
fn = inspect.getfile(get_this_dir)
current_dir = os.path.dirname(fn)

load_file = os.path.join(current_dir, "begin.py")
load_file = os.path.normpath(load_file)
if (os.path.isfile(load_file)):
    exec(compile(open(load_file, "rb").read(), load_file, 'exec'), globals())

log = StreamIO(sys.stderr, sys.stdout)

result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)

if (unittest_output):
    print("\n")
    print(unittest_output.getvalue())
    unittest_output.close()
else:
    print("BL ERROR:: no unittest output")

if (have_test_skip):
    print("\nUnittest skip exists!")
    if (result.skipped):
        print("...with the following skipping message(s):")
        for skipped in result.skipped:
            print("   ", skipped[1])
        print("\n")
    else:
        print("No tests skipped")
else:
    print("\nUnitest skip does not exist!")
    if (skipped_tests):
        print("The following tests were skipped and marked as passed:")
        for test in skipped_tests:
            print("   ", test)
        print("\n")

# finish off closing all maps (maybe use end.py file?!)
list(map(coot.close_molecule, coot_utils.molecule_number_list()))

# some garbage collection. Not sure if this will really shed light into
# reference counting e.g.
print("\nFINALLY collect some GARBAGE\n")
import gc
import inspect
gc.enable()
gc.set_debug(gc.DEBUG_LEAK)
print("BL INFO:: no of collected carbage items %s\n" %gc.collect())
all_garbage = gc.garbage
print("BL DEBUG:: list of garbage items:\n")
for garb in all_garbage:
    print("    ", garb)
print("\nDetailed information for functions:\n")
for garb in all_garbage:
    if inspect.isfunction(garb) and not inspect.isbuiltin(garb):
        print("  BL INFO:: garbage item and location:")
        print("    %s" %garb)
        print("    %s" %gc.get_referents(garb)[0])
        print()

# cheating?! We only exit Coot if we are not in graphics mode
if (coot.use_graphics_interface_state() == 0):
    if (result.wasSuccessful()):
        coot.coot_real_exit(0)
    else:
        coot.coot_real_exit(1)
