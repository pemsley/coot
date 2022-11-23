# 02_shelx.py
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

import coot
import coot_utils
import begin
import unittest
import os

# NOTE: hof.fcf has  _refln_A_calc and _refln_B_calc, not fcalc and phase.
# This will not work.  Why is it here?
#
hof_fcf = os.path.join(begin.unittest_data_dir, "hof.fcf")
hof_res = os.path.join(begin.unittest_data_dir, "HOF.RES")

#global imol_hof_res
imol_hof_res = False

insulin_fcf = os.path.join(begin.unittest_data_dir, "insulin.fcf")
insulin_res = os.path.join(begin.unittest_data_dir, "insulin.res")
hollander_ins = os.path.join(begin.unittest_data_dir, "hollander.ins")
global imol_insulin_res
global imol_insulin_map
imol_insulin_res = -1 # see later
imol_insulin_map = -1 # see later
m_miller_res = os.path.join(begin.unittest_data_dir, "miller/shelx-test4-NPD-mini.res")


class ShelxTestFunctions(unittest.TestCase):

    def test01_0(self):
        """Read small molecule .res file"""

        imol = coot.read_pdb(hof_res)
        self.assertTrue(coot_utils.valid_model_molecule_qm(imol))
        imol_hof_res = imol


    def test02_0(self):
        """Read hollander small molecule .res file"""

        if self.skip_test(not os.path.isfile(hollander_ins),
                          "%s  does not exist - skipping test" %hollander_ins):
            return

        imol = coot.read_pdb(hollander_ins)
        self.assertTrue(coot_utils.valid_model_molecule_qm(imol), "   fail: bad molecule for %s" %hollander_ins)
        spg = coot.show_spacegroup(imol)
        self.assertEqual(spg, "I 41 2 2", "   fail: wrong spacegroup for %s %s" %(hollander_ins, spg))


    def test03_0(self):
        """read shelx insulin with fcf"""

        global imol_insulin_res
        global imol_insulin_map

        imol_insulin_res_local = coot.handle_read_draw_molecule_with_recentre(insulin_res(), 1)
        self.assertTrue(coot_utils.valid_model_molecule_qm(imol_insulin_res_local),
                        "   Bad insulin.res: %s for %s" %(insulin_res, imol_insulin_res_local))

        imol_insulin_res = imol_insulin_res_local  # used in water addition test
        imol = coot.handle_shelx_fcf_file(insulin_fcf)
        self.assertTrue(coot_utils.valid_map_molecule_qm(imol), "    Bad read of %s %s" %(insulin_fcf, imol))

        # remove the name test for Sigmaa now that we use read_small_molecule_data_cif

        self.assertEqual(coot.show_spacegroup(imol),
                             coot.show_spacegroup(imol_insulin_res_local),
                             "   Mismatch spacegroups %s %s" %(coot.show_spacegroup(imol),
                                                               coot.show_spacegroup(imol_insulin_res_local)))

        self.assertEqual(coot.show_spacegroup(imol), "I 21 3",
                             "   Bad spacegroups %s" %(coot.show_spacegroup(imol)))

        # good then
        imol_insulin_map = imol
        coot.rotate_y_scene(coot.rotate_n_frames(200), 0.1)
        # close_molecule(imol)    # needed later
        # imol_insulin_res_local needed later -> global?

    # The "SHELX Woe" problem
    #
    def test04_0(self):
        """Write an INS from PDB test"""

        # First, check return status on a bogus molecule
        status = coot.write_shelx_ins_file(204050, "no_molecule.ins")
        self.assertEqual(status, 0,
                             "bad exit status from write_shelx_ins_file on bogus molecule %s" %status)

        # Now check return status on a bogus molecule.  The Happy Path.
        #
        self.assertTrue(coot_utils.valid_model_molecule_qm(begin.imol_rnase),
                        "   imol_rnase not valid.")
        rnase_ins = "rnase.ins"
        status = coot.write_shelx_ins_file(begin.imol_rnase, rnase_ins)
        self.assertEqual(status, 1, "   failure to write INS file %s from PDB: status %s" %(rnase_ins, status))


    def test05_0(self):
        """new molecule by atom selection inherits shelx molecule flag"""

        global imol_insulin_res
        insulin_frag = coot.new_molecule_by_atom_selection(imol_insulin_res, "//B/2010-2020")
        self.assertTrue(coot_utils.valid_model_molecule_qm(insulin_frag),
                        " bad fragment of insulin res molecule")
        self.assertTrue(coot_utils.shelx_molecule_qm(insulin_frag),
                        " bad shelx flag from insulin-frag")


    def test06_0(self):
        """Addition of Terminal Residue on SHELX molecule has correct occupancy"""

        global imol_insulin_res
        global imol_insulin_map
        insulin_frag = coot.new_molecule_by_atom_selection(imol_insulin_res, "//B/2010-2020")
        self.assertTrue(coot_utils.valid_model_molecule_qm(insulin_frag),
                        " bad fragment of insulin res molecule")

        coot.set_imol_refinement_map(imol_insulin_map)
        coot.add_terminal_residue(insulin_frag, "B", 2020, "ALA", 1)

        res_atoms = coot.residue_info(insulin_frag, "B", 2021, "")

        self.assertTrue(res_atoms,
                        " bad residue info after add terminal residue")

        for atom in res_atoms:
            occ = atom[1][0]
            self.assertAlmostEqual(occ, 11.0, 1, " bad occupancides in new residue %s" %atom)


    def test07_0(self):
        """Add water to SHELX molecule"""

        global imol_insulin_res
        self.assertTrue(coot_utils.valid_model_molecule_qm(imol_insulin_res),
                        "   failed to get a valid imol_insulin_res")

        coot.set_pointer_atom_molecule(imol_insulin_res)
        coot.set_rotation_centre(3, -1, 60)
        coot.place_typed_atom_at_pointer("Water")
        # test is to have to occupancy of the new HOH to be 11.0
        coot_utils.shelx_waters_all_good_occ_qm(self, imol_insulin_res)


    # Tobias Beck test
    def test08_0(self):
        """Find Waters for a SHELXL molecule"""

        self.assertTrue(coot_utils.valid_model_molecule_qm(imol_insulin_res),
                        "   failed to get a valid imol_insulin_res")

        n_chains_pre = coot.n_chains(imol_insulin_res)
        coot.find_waters(imol_insulin_map, imol_insulin_res, 0, 0.6, 1)
        n_chains_post = coot.n_chains(imol_insulin_res)
        self.assertTrue(n_chains_pre == n_chains_post,
                        "Find waters on a shelx molecule created a new chain %s %s"
                        %(n_chains_pre, n_chains_post))
        coot_utils.shelx_waters_all_good_occ_qm(self, imol_insulin_res)


# non positive definite anistropic atom (reported by Mitch Miller)
# crash test
    def test09_0(self):
        """NPD Anisotropic Atom [Mitch Miller]"""

        imol_miller = coot.handle_read_draw_molecule_with_recentre(m_miller_res, 1)

        self.assertTrue(coot_utils.valid_model_molecule_qm(imol_miller), "Bad read of miller test molecule")

        coot.set_show_aniso(1)  # crash?
        coot.rotate_y_scene(coot.rotate_n_frames(100), 0.1)
        coot.close_molecule(imol_miller)
        coot.set_show_aniso(0)

# cheesy test, close the shelx molecules

    def test10_0(self):
        """close shelx molecules"""

        coot.close_molecule(imol_insulin_map)
        coot.close_molecule(imol_insulin_res)

        self.assertFalse(coot_utils.valid_model_molecule_qm(imol_insulin_res),
                    "imol_insulin_res: %s is still valid after closure!"
                    %imol_insulin_res)
        self.assertFalse(coot_utils.valid_map_molecule_qm(imol_insulin_map),
                    "imol_insulin_map: %s is still valid after closure!"
                    %imol_insulin_map)


    def test11_0(self):
        """Aniso Bs in P21"""

        from types import ListType

        def aniso_b_from_atom(atom):
            occ_etc = atom[1]
            return occ_etc[1]

        imol = begin.unittest_pdb("horma-p21.res")
        self.assertTrue(coot_utils.valid_model_molecule_qm(imol),
                        "  failed to get a valid imoll from horma-p21.res")
        coot.write_shelx_ins_file(imol, "new-horma.ins")
        imol_2 = coot.read_pdb("new-horma.ins")
        at_1 = coot_utils.get_atom(imol   , "A", 4, "", " N1 ")
        at_2 = coot_utils.get_atom(imol_2 , "A", 4, "", " N1 ")
        b_1  = aniso_b_from_atom(at_1)
        b_2  = aniso_b_from_atom(at_2)

        print("b_1:", b_1)
        print("b_2:", b_2)

        self.assertTrue(type(b_1) is ListType)
        self.assertTrue(type(b_2) is ListType)

        list(map(lambda x, y: self.assertAlmostEqual(x, y), b_1, b_2))


    def test12_0(self):
        """Don't crash on reading a strange HAT file"""

        coot.handle_read_draw_molecule_with_recentre(os.path.join(begin.unittest_data_dir(), "crash.hat"), 0)
        # this is nonsense since we already 'tested' for read/crash
        # but it is a proper test...
        imol = begin.unittest_pdb("crash.hat")
        self.assertTrue(coot_utils.valid_model_molecule_qm(imol))
