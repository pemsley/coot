# unittest_extensions.py
# Copyright 2007, 2008, 2009 by The University of York
# Copyright 2008, 2009 by Bernhard Lohkamp
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

import os
home_dir = os.getenv('HOME')
if (not home_dir and os.name == 'nt'):
        home_dir = os.getenv('COOT_HOME')

project_dir = os.path.join(home_dir, "Projects", "coot", "python-tests")
unittest_filename = os.path.join(project_dir, "coot_unittest.py")
begin_filename = os.path.join(project_dir, "begin.py")


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
            # search, starting from HOME
            if (os.name == 'nt'):
                home = os.getenv("COOT_HOME")
            else:
                home = os.getenv("HOME")
            # now we assume its in ~/Projects/coot
            full_dir = os.path.normpath(os.path.join(home, "Projects", "coot", unittest_dir))
            if (os.path.isdir(full_dir)):
                full_name = os.path.join(full_dir, unittest_file)
                if (os.path.isfile(full_name)):
                    execfile(full_name, globals())
                else:
                    print "BL INFO:: could not find %s in %s" %(unittest_file, full_dir)
            else:
                print "BL INFO:: searching $HOME..."
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

        menu = coot_gui.coot_menubar_menu("Unittesting")


        def run_all_tests_func():
            exec_file("coot_unittest.py")

        coot_gui.add_simple_coot_menu_menuitem(menu, "Run All Unittest Tests",
                                      lambda func: run_all_tests_func())


        def run_test_set_gui():
            exec_file("begin.py")
            #print "BL DEBUG:: test list", test_list
            #for tmp in test_list:
            #    print "BL DEBUG:: test str", str(tmp).split(".")[1].rstrip("\'>")

            # make buttons list:
            buttons = []
            for test in test_list:
                buttons.append([str(test).split(".")[1].rstrip("\'>"), "run_test_set(" + str(test_list.index(test)) + ")"])

            #run_test_set(7)

            coot_gui.dialog_box_of_buttons("Coot Unittests Sets", [200, 320], buttons,
                                  "  Close  ")

            

        coot_gui.add_simple_coot_menu_menuitem(menu, "Run A Unittest Set...",
                                      lambda func: run_test_set_gui())

        def run_one_test_gui():
            exec_file("begin.py")
            test_list = list_of_all_tests()

            # make buttons list:
            buttons = []
            prev_test_class = "PdbMtzTestFunction"
            for test in test_list:
              if (not prev_test_class in test[2]):
                # new set, let's introduce
                #print "BL DEBUG:: test id of new set", test[2]
                buttons.append(["HSep"])
              buttons.append([test[1], "run_one_test(" + str(test_list.index(test)) + ")"])
              prev_test_class = test[2].split(".")[1]
            
            coot_gui.dialog_box_of_buttons("Coot Unittests", [400, 500], buttons,
                                  "  Close  ")
            

        coot_gui.add_simple_coot_menu_menuitem(menu, "Run One Unittest Test...",
                                      lambda func: run_one_test_gui())

