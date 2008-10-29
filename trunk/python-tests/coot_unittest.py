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

def get_this_dir():
    pass


# get directory of this file and execute tests found in this dir
fn = inspect.getfile(get_this_dir)
current_dir = os.path.dirname(fn)

load_file = os.path.join(current_dir, "begin.py")
load_file = os.path.normpath(load_file)
if (os.path.isfile(load_file)):
    execfile(load_file, globals())

log = StreamIO(sys.stderr, sys.stdout)

result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)

if (unittest_output):
    print "\n"
    print unittest_output.getvalue()
    unittest_output.close()
else:
    print "BL ERROR:: no unittest output"

if (have_test_skip):
    print "\nUnittest skip exists!"
    if (result.skipped):
        print "...with the following skipping message(s):"
        for skipped in result.skipped:
            print "   ", skipped[1]
        print "\n"
    else:
        print "No tests skipped"
else:
    print "\nUnitest skip does not exist!"
    if (skipped_tests):
        print "The following tests were skipped and marked as passed:"
        for test in skipped_tests:
            print "   ", test
        print "\n"

# cheating?! We only exit Coot if we are not in graphics mode
if (use_graphics_interface_state() == 0):
    if (result.wasSuccessful()):
        coot_real_exit(0)
    else:
        coot_real_exit(1)

