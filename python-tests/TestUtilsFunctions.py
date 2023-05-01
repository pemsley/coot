# 08_utils.py
# Copyright 2008 by The University of York
# Author: Bernhard Lohkamp
# Copyright 2007, 2008 by The University of Oxford
# Author: Paul Emsley

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

import unittest
import os
import coot
import coot_utils
import coot_testing_utils


class TestUtilsFunctions(unittest.TestCase):

    def test01_0(self):
        """Test key symbols"""

        # coot.add_key_binding_gtk3_py("name", "missing", "") # old
        # coot.add_key_binding_gtk3_py(key-code, Ctrl-mode, function, label)

        test_list = [[coot.key_sym_code_py("a-symbol"), -1],
                     [coot.key_sym_code_py(["a", "b"]), -1],
                     [coot.key_sym_code_py("A"), 65],
                     [coot.key_sym_code_py("a"), 97],
                     [coot.key_sym_code_py("9"), 57],
                     [coot.key_sym_code_py("cent"), 162],
                     [coot.key_sym_code_py(">"), 62],
                     [coot.key_sym_code_py(";"), 59],
                     [coot.key_sym_code_py("|"), 124],
                     ]

        for code, result in test_list:
            self.assertEqual(code, result, " fail on key_sym_code, %s is not equal to %s " %(code, result))

    # 20230213-PE No scheme today sadly
    # def test02_0(self):
    #     """Test running a scheme function"""

    #     tot = coot.run_scheme_command("(+ 2 4)")
    #     self.assertEqual(tot, 6)

    #     #run_scheme_command("define test-val 4")
    #     #run_scheme_command("(set! test-val 2)")
    #     rv = coot.run_scheme_command("(rotation-centre)")
    #     print("BL DEBUG:: return scheme is ", rv)
    #     #rv = run_scheme_command("test-val")
    #     #print "BL DEBUG:: return scheme is ", rv
    #     rv = coot.run_scheme_command("2")
    #     self.assertEqual(rv, 2)


    def test03_0(self):
        """Internal/External Molecule Numbers match"""

        m = coot_utils.molecule_number_list()
        #print "   m: ", m
        #print " own: ", map(own_molecule_number, m)
        self.assertTrue(m == list(map(coot.own_molecule_number, m)))


        def test04_0(self):
            """spacegroup operators to space group conversion"""

            test_groups = ["P 1", "P -1", "P 1 2 1", "P 1 1 2", "P 2 1 1", "P 1 21 1", "P 1 1 21",
                           "P 21 1 1", "C 1 2 1", "A 1 2 1", "I 1 2 1", "A 1 1 2", "B 1 1 2", "I 1 1 2",
                           "B 2 1 1", "C 2 1 1", "I 2 1 1", "P 1 m 1", "P 1 1 m", "P m 1 1", "P 1 c 1",
                           "P 1 n 1", "P 1 a 1", "P 1 1 a", "P 1 1 n", "P 1 1 b", "P b 1 1", "P n 1 1",
                           "P c 1 1", "C 1 m 1", "A 1 m 1", "I 1 m 1", "A 1 1 m", "B 1 1 m", "I 1 1 m",
                           "B m 1 1", "C m 1 1", "I m 1 1", "C 1 c 1", "A 1 n 1", "I 1 a 1", "A 1 a 1",
                           "C 1 n 1", "I 1 c 1", "A 1 1 a", "B 1 1 n", "I 1 1 b", "B 1 1 b", "A 1 1 n", "I 1 1 a",
                           "B b 1 1", "C n 1 1", "I c 1 1", "C c 1 1", "B n 1 1", "I b 1 1", "P 1 2/m 1", "P 1 1 2/m",
                           "P 2/m 1 1", "P 1 21/m 1", "P 1 1 21/m", "P 21/m 1 1", "C 1 2/m 1", "A 1 2/m 1", "I 1 2/m 1",
                           "A 1 1 2/m", "B 1 1 2/m", "I 1 1 2/m", "B 2/m 1 1", "C 2/m 1 1", "I 2/m 1 1", "P 1 2/c 1",
                           "P 1 2/n 1", "P 1 2/a 1", "P 1 1 2/a", "P 1 1 2/n", "P 1 1 2/b", "P 2/b 1 1", "P 2/n 1 1",
                           "P 2/c 1 1", "P 1 21/c 1", "P 1 21/n 1", "P 1 21/a 1", "P 1 1 21/a", "P 1 1 21/n", "P 1 1 21/b",
                           "P 21/b 1 1", "P 21/n 1 1", "P 21/c 1 1", "C 1 2/c 1", "A 1 2/n 1", "I 1 2/a 1", "A 1 2/a 1",
                           "C 1 2/n 1", "I 1 2/c 1", "A 1 1 2/a", "B 1 1 2/n", "I 1 1 2/b", "B 1 1 2/b", "A 1 1 2/n",
                           "I 1 1 2/a", "B 2/b 1 1", "C 2/n 1 1", "I 2/c 1 1", "C 2/c 1 1", "B 2/n 1 1", "I 2/b 1 1",
                           "P 2 2 2", "P 2 2 21", "P 21 2 2", "P 2 21 2", "P 21 21 2", "P 2 21 21", "P 21 2 21",
                           "P 21 21 21", "C 2 2 21", "A 21 2 2", "B 2 21 2", "C 2 2 2", "A 2 2 2", "B 2 2 2",
                           "F 2 2 2", "I 2 2 2", "I 21 21 21", "P m m 2", "P 2 m m", "P m 2 m", "P m c 21", "P c m 21",
                           "P 21 m a", "P 21 a m", "P b 21 m", "P m 21 b", "P c c 2", "P 2 a a", "P b 2 b", "P m a 2",
                           "P b m 2", "P 2 m b", "P 2 c m", "P c 2 m", "P m 2 a", "P c a 21", "P b c 21", "P 21 a b",
                           "P 21 c a", "P c 21 b", "P b 21 a", "P n c 2", "P c n 2", "P 2 n a", "P 2 a n", "P b 2 n",
                           "P n 2 b", "P m n 21", "P n m 21", "P 21 m n", "P 21 n m", "P n 21 m", "P m 21 n",
                           "P b a 2", "P 2 c b", "P c 2 a", "P n a 21", "P b n 21", "P 21 n b", "P 21 c n",
                           "P c 21 n", "P n 21 a", "P n n 2", "P 2 n n", "P n 2 n", "C m m 2", "A 2 m m", "B m 2 m",
                           "C m c 21", "C c m 21", "A 21 m a", "A 21 a m", "B b 21 m", "B m 21 b", "C c c 2", "A 2 a a",
                           "B b 2 b", "A m m 2", "B m m 2", "B 2 m m", "C 2 m m", "C m 2 m", "A m 2 m", "A b m 2",
                           "B m a 2", "B 2 c m", "C 2 m b", "C m 2 a", "A c 2 m", "A m a 2", "B b m 2", "B 2 m b", "C 2 c m",
                           "C c 2 m", "A m 2 a", "A b a 2", "B b a 2", "B 2 c b", "C 2 c b", "C c 2 a", "A c 2 a", "F m m 2",
                           "F 2 m m", "F m 2 m", "F d d 2", "F 2 d d", "F d 2 d", "I m m 2", "I 2 m m", "I m 2 m", "I b a 2",
                           "I 2 c b", "I c 2 a", "I m a 2", "I b m 2", "I 2 m b", "I 2 c m", "I c 2 m", "I m 2 a",
                           "P m m m"   "P c c m", "P m a a", "P b m b",
                           "P m m a", "P m m b",
                           "P b m m", "P c m m", "P m c m", "P m a m", "P n n a", "P n n b", "P b n n", "P c n n",
                           "P n c n", "P n a n", "P m n a", "P n m b", "P b m n", "P c n m", "P n c m", "P m a n",
                           "P c c a", "P c c b", "P b a a", "P c a a", "P b c b", "P b a b", "P b a m", "P m c b",
                           "P c m a", "P c c n", "P n a a", "P b n b", "P b c m", "P c a m", "P m c a", "P m a b",
                           "P b m a", "P c m b", "P n n m", "P m n n", "P n m n",
                           "P b c n", "P c a n", "P n c a", "P n a b", "P b n a",
                           "P c n b", "P b c a", "P c a b", "P n m a", "P m n b", "P b n m", "P c m n", "P m c n",
                           "P n a m", "C m c m", "C c m m", "A m m a", "A m a m", "B b m m", "B m m b", "C m c a",
                           "C c m b", "A b m a", "A c a m", "B b c m", "B m a b", "C m m m", "A m m m", "B m m m",
                           "C c c m", "A m a a", "B b m b", "C m m a", "C m m b", "A b m m", "A c m m", "B m c m",
                           "B m a m",
                           "F m m m"  "I m m m", "I b a m", "I m c b",
                           "I c m a", "I b c a", "I c a b", "I m m a", "I m m b", "I b m m", "I c m m", "I m c m",
                           "I m a m", "P 4", "P 41", "P 42", "P 43", "I 4", "I 41", "P -4", "I -4", "P 4/m", "P 42/m",
                           "I 4/m",
                           "P 4 2 2", "P 4 21 2", "P 41 2 2", "P 41 21 2", "P 42 2 2", "P 42 21 2", "P 43 2 2",
                           "P 43 21 2", "I 4 2 2", "I 41 2 2", "P 4 m m", "P 4 b m", "P 42 c m", "P 42 n m",
                           "P 4 c c", "P 4 n c", "P 42 m c", "P 42 b c", "I 4 m m", "I 4 c m", "I 41 m d", "I 41 c d",
                           "P -4 2 m", "P -4 2 c", "P -4 21 m", "P -4 21 c", "P -4 m 2", "P -4 c 2", "P -4 b 2",
                           "P -4 n 2", "I -4 m 2", "I -4 c 2", "I -4 2 m", "I -4 2 d", "P 4/m m m", "P 4/m c c",
                           "P 4/m b m", "P 4/m n c",
                           "P 42/m m c",
                           "P 42/m c m",
                           "P 42/m b c", "P 42/m n m",
                           "I 4/m m m", "I 4/m c m",
                           "P 3", "P 31", "P 32", "P -3",
                           "P 3 1 2", "P 3 2 1", "P 31 1 2", "P 31 2 1", "P 32 1 2", "P 32 2 1",
                           "P 3 m 1", "P 3 1 m", "P 3 c 1", "P 3 1 c",
                           "P -3 1 m", "P -3 1 c", "P -3 m 1", "P -3 c 1",
                           "P 6", "P 61", "P 65", "P 62", "P 64", "P 63", "P -6", "P 6/m",
                           "P 63/m", "P 6 2 2", "P 61 2 2", "P 65 2 2", "P 62 2 2", "P 64 2 2", "P 63 2 2",
                           "P 6 m m", "P 6 c c", "P 63 c m", "P 63 m c", "P -6 m 2", "P -6 c 2", "P -6 2 m",
                           "P -6 2 c", "P 6/m m m", "P 6/m c c", "P 63/m c m", "P 63/m m c", "P 2 3", "F 2 3",
                           "I 2 3", "P 21 3", "I 21 3", "P m -3", "F m -3",
                           "I m -3", "P a -3", "I a -3", "P 4 3 2", "P 42 3 2", "F 4 3 2", "F 41 3 2",
                           "I 4 3 2", "P 43 3 2", "P 41 3 2", "I 41 3 2", "P -4 3 m", "F -4 3 m", "I -4 3 m", "P -4 3 n",
                           "F -4 3 c", "I -4 3 d", "P m -3 m", "P m -3 n",
                           "F m -3 m", "F m -3 c",
                           "I m -3 m", "I a -3 d"]

            # failures: "P n n n :1" "P n n n :2"  "P b a n :1"  "P b a n :2" 
            # "P n c b :1" "P n c b :2" "P c n a :1" "P c n a :2"  "P m m n :1"
            # "P m m n :2"  "P n m m :1"  "P n m m :2" "P m n m :1" "P m n m :2" 
            # "R 3 :H"  "R 3 :R" "R -3 :H" "R -3 :R" "R 3 2 :H"  "R 3 2 :R" 
            # 
            # probable failures (ran out of patience to test the indiviually)
            # "C c c a :1" "C c c a :2" "C c c b :1" "C c c b :2" "A b a a :1"
            # "A b a a :2" "A c a a :1" "A c a a :2" "B b c b :1" "B b c b :2" "B b a b :1"
            # "B b a b :2" "F d d d :1" "F d d d :2"
            # "P 4/n :1" "P 4/n :2" "P 42/n :1" "P 42/n :2"
            # "I 41/a :1" "I 41/a :2" 
            # "P 4/n b m :1" "P 4/n b m :2" "P 4/n n c :1" "P 4/n n c :2" 
            #  "P 4/n m m :1" "P 4/n m m :2" "P 4/n c c :1" "P 4/n c c :2" 
            #  "P 42/n b c :1" "P 42/n b c :2" "P 42/n n m :1" "P 42/n n m :2" 
            # "P 42/n m c :1" "P 42/n m c :2" "P 42/n c m :1" "P 42/n c m :2" 
            # "I 41/a m d :1" "I 41/a m d :2" "I 41/a c d :1" "I 41/a c d :2"
            #  "P n -3 :1" "P n -3 :2"
            # "F d -3 :1"
            # "F d -3 m :1" "F d -3 m :2" "F d -3 c :1"
            # "P n -3 m :1" "F d -3 :2"
            # "R 3 m :H" "R 3 m :R" "R 3 c :H"  "R 3 c :R"  "R -3 m :H" "R -3 m :R"
            #  "R -3 c :H" "R -3 c :R" "P n -3 n :1" "P n -3 n :2" 
            # "P n -3 m :2"
            # "F d -3 c :2" 

            # clipper doesn't know about these space groups (mmdb returns
            # the symop strings correctly (below called symops) (in that
            # they correspond to syminfo.lib):
            #
            # "I 1 21 1" "C 1 21 1"  "B 1 1 m"

            from types import ListType
            from types import StringType

            imol = coot_testing_utils.unittest_pdb("monomer-ACT.pdb")

            for space_group in test_groups:

                set_success = coot.set_space_group(imol, space_group)
                symops = coot.symmetry_operators(imol)

                self.assertTrue(set_success == 1,
                                "   bad status on setting space group %s" %space_group)

                self.assertTrue(type(symops) is ListType,
                                "   bad symops for %s: %s" %(space_group, symops))

                derived_HM = coot.symmetry_operators_to_xHM(symops)

                self.assertTrue(type(derived_HM) is StringType,
                                "   bad derived HM %s" %derived_HM)

                self.assertTrue(space_group == derived_HM,
                                "   No match %s and %s" %(space_group, derived_HM))

