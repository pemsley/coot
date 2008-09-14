# unittest_extensions.py
# Copyright 2007, 2008 by The University of York
# Copyright 2008 by Bernhard Lohkamp
# Copyright 2007 by Paul Emsley
# Copyright 2007 by The University of Oxford
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

unittest_filename="/y/people/lohkamp/Projects/coot/python-tests/coot_unittest.py"
begin_filename="/y/people/lohkamp/Projects/coot/python-tests/begin.py"

have_coot_python = False
try: 
  import coot_python
  have_coot_python = True
except:
    pass

def exec_file(filename):
    if (os.path.isfile(unittest_filename) and
        (filename == "coot_unittest.py")):
        execfile(unittest_filename, globals())
    elif (os.path.isfile(begin_filename) and
          (filename == "begin.py")):
        execfile(begin_filename, globals())
    else:
        unittest_dir = "python-tests"
        unittest_file = filename
        # assume we started coot in bin:
        full_dir = os.path.normpath(os.path.join("..", unittest_dir))
        if (os.path.isdir(full_dir)):
            full_name = os.path.join(unittest_dir, unittest_file)
            if (os.path.isfile(full_name)):
                execfile(full_name, globals())
            else:
                print "BL INFO:: could not find %s in ../%s" %(unittest_file, unittest_dir)
        else:
            print "BL INFO:: could not find directory ../%s" %unittest_dir
            print "BL INFO:: searching $HOME..."
            # search, starting from HOME
            if (os.name == 'nt'):
                home = os.getenv("COOT_HOME")
            else:
                home = os.getenv("HOME")
            for root, dirs, files in os.walk(home):
                #print "BL DEBUG:: root, dirs", root, dirs
                if unittest_dir in dirs:
                    full_name = os.path.join(root, unittest_dir, unittest_file)
                    if (os.path.isfile(full_name)):
                        execfile(full_name, globals())
                        break
                    else:
                        print "BL INFO:: could not find %s in %s/%s" %(unittest_file, root, unittest_dir)
                        print "BL INFO:: continue searching..."

if (have_coot_python):
    if coot_python.main_menubar():

        menu = coot_menubar_menu("Unittesting")


        def run_all_tests_func():
            exec_file("coot_unittest.py")

        add_simple_coot_menu_menuitem(menu, "Run All Unittest Tests",
                                      lambda func: run_all_tests_func())


        def run_test_set_gui():
            exec_file("begin.py")
            print "BL DEBUG:: test list", test_list
            for tmp in test_list:
                print "BL DEBUG:: test str", str(tmp).split(".")[1].rstrip("\'>")

            # make buttons list:
            buttons = []
            for test in test_list:
                buttons.append([str(test).split(".")[1].rstrip("\'>"), "run_test_set(" + str(test_list.index(test)) + ")"])

            #run_test_set(7)

            dialog_box_of_buttons("Coot Unittests Sets", [200, 320], buttons,
                                  "  Close  ")

            

        add_simple_coot_menu_menuitem(menu, "Run A Unittest Set...",
                                      lambda func: run_test_set_gui())

        def run_one_test_gui():
            exec_file("begin.py")
            test_list = list_of_all_tests()

            # make buttons list:
            buttons = []
            for test in test_list:
                buttons.append([test[1], "run_one_test(" + str(test_list.index(test)) + ")"])
            
            dialog_box_of_buttons("Coot Unittests", [400, 500], buttons,
                                  "  Close  ")
            

        add_simple_coot_menu_menuitem(menu, "Run One Unittest Test...",
                                      lambda func: run_one_test_gui())

