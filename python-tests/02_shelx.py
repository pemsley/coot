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

import unittest
import os

# NOTE: hof.fcf has  _refln_A_calc and _refln_B_calc, not fcalc and phase.
# This will not work.  Why is it here?
#
hof_fcf = os.path.join(unittest_data_dir, "hof.fcf")
hof_res = os.path.join(unittest_data_dir, "HOF.RES")

#global imol_hof_res
imol_hof_res = False

insulin_fcf = os.path.join(unittest_data_dir, "insulin.fcf")
insulin_res = os.path.join(unittest_data_dir, "insulin.res")
hollander_ins = os.path.join(unittest_data_dir, "hollander.ins")
global imol_insulin_res
global imol_insulin_map
imol_insulin_res = -1 # see later
imol_insulin_map = -1 # see later
m_miller_res = os.path.join(unittest_data_dir, "miller/shelx-test4-NPD-mini.res")


class ShelxTestFunctions(unittest.TestCase):

    def test01_0(self):
        """Read small molecule .res file"""

        if self.skip_test(not type(hof_res) is StringType,
                          "hof-res not defined - skipping test"):
            return
        if self.skip_test(not os.path.isfile(hof_res),
                          "%s does not exist - skipping test" %hof_res):
            return
            
        imol = read_pdb(hof_res)
        self.failUnless(valid_model_molecule_qm(imol))
        imol_hof_res = imol


    def test02_0(self):
        """Read hollander small molecule .res file"""

        if self.skip_test(not os.path.isfile(hollander_ins),
                          "%s  does not exist - skipping test" %hollander_ins):
            return

        imol = read_pdb(hollander_ins)
        self.failUnless(valid_model_molecule_qm(imol),
                        "   fail: bad molecule for %s" %hollander_ins)
        spg = show_spacegroup(imol)
        self.failUnlessEqual(spg, "I 41 2 2",
                             "   fail: wrong spacegroup for %s %s" %(hollander_ins, spg))
        

    def test03_0(self):
        """read shelx insulin with fcf"""

        global imol_insulin_res
        global imol_insulin_map
        
        imol_insulin_res_local = handle_read_draw_molecule_with_recentre(insulin_res(), 1)
        self.failUnless(valid_model_molecule_qm(imol_insulin_res_local),
                        "   Bad insulin.res: %s for %s" %(insulin_res, imol_insulin_res_local))

        imol_insulin_res = imol_insulin_res_local  # used in water addition test
        imol = handle_shelx_fcf_file(insulin_fcf)
        self.failUnless(valid_map_molecule_qm(imol),
                        "    Bad read of %s %s" %(insulin_fcf, imol))
        name = molecule_name(imol)
        cif_name = insulin_fcf + ".cif SigmaA"
        self.failUnlessEqual(name, cif_name,
                             "   Bad name match %s != %s" %(name, cif_name))

        self.failUnlessEqual(show_spacegroup(imol), show_spacegroup(imol_insulin_res_local),
                             "   Mismatch spacegroups %s %s" %(show_spacegroup(imol),
                                                               show_spacegroup(imol_insulin_res_local)))

        self.failUnlessEqual(show_spacegroup(imol), "I 21 3",
                             "   Bad spacegroups %s" %(show_spacegroup(imol)))

        # good then
        imol_insulin_map = imol
        rotate_y_scene(rotate_n_frames(200), 0.1)
        # close_molecule(imol)    # needed later
        # imol_insulin_res_local needed later -> global?
        

    # The "SHELX Woe" problem
    #
    def test04_0(self):
        """Write an INS from PDB test"""

        # First, check return status on a bogus molecule
        status = write_shelx_ins_file(204050, "no_molecule.ins")
        self.failUnlessEqual(status, 0,
                             "bad exit status from write_shelx_ins_file on bogus molecule %s" %status)

        # Now check return status on a bogus molecule.  The Happy Path.
        #
        self.failUnless(valid_model_molecule_qm(imol_rnase),
                        "   imol-rnase not valid.")
        rnase_ins = "rnase.ins"
        status = write_shelx_ins_file(imol_rnase, rnase_ins)
        self.failUnlessEqual(status, 1,
                             "   failure to write INS file %s from PDB: status %s" %(rnase_ins, status))


    def test05_0(self):
        """new molecule by atom selection inherits shelx molecule flag"""

        global imol_insulin_res
        insulin_frag = new_molecule_by_atom_selection(imol_insulin_res,
                                                      "//B/2010-2020")
        self.failUnless(valid_model_molecule_qm(insulin_frag),
                        " bad fragment of insulin res molecule")
        self.failUnless(shelx_molecule_qm(insulin_frag),
                        " bad shelx flag from insulin-frag")


    def test06_0(self):
        """Addition of Terminal Residue on SHELX molecule has correct occupancy"""

        global imol_insulin_res
        global imol_insulin_map
        insulin_frag = new_molecule_by_atom_selection(imol_insulin_res,
                                                      "//B/2010-2020")
        self.failUnless(valid_model_molecule_qm(insulin_frag),
                        " bad fragment of insulin res molecule")

        set_imol_refinement_map(imol_insulin_map)
        add_terminal_residue(insulin_frag, "B", 2020, "ALA", 1)

        res_atoms = residue_info(insulin_frag, "B", 2021, "")

        self.failUnless(res_atoms,
                        " bad residue info after add terminal residue")

        for atom in res_atoms:
            occ = atom[1][0]
            self.failUnlessAlmostEqual(occ, 11.0, 1,
                                       " bad occupancides in new residue %s" %atom)
        

    def test07_0(self):
        """Add water to SHELX molecule"""

        global imol_insulin_res
        self.failUnless(valid_model_molecule_qm(imol_insulin_res),
                        "   failed to get a valid imol_insulin_res")
        
        set_pointer_atom_molecule(imol_insulin_res)
        set_rotation_centre(3, -1, 60)
        place_typed_atom_at_pointer("Water")
        # test is to have to occupancy of the new HOH to be 11.0
        shelx_waters_all_good_occ_qm(self, imol_insulin_res)


    # Tobias Beck test
    def test08_0(self):
        """Find Waters for a SHELXL molecule"""

        self.failUnless(valid_model_molecule_qm(imol_insulin_res),
                        "   failed to get a valid imol_insulin_res")

        n_chains_pre = n_chains(imol_insulin_res)
        find_waters(imol_insulin_map, imol_insulin_res, 0, 0.6, 1)
        n_chains_post = n_chains(imol_insulin_res)
        self.failUnless(n_chains_pre == n_chains_post,
                        "Find waters on a shelx molecule created a new chain %s %s"
                        %(n_chains_pre, n_chains_post))
        shelx_waters_all_good_occ_qm(self, imol_insulin_res)
        

# non positive definite anistropic atom (reported by Mitch Miller)
# crash test
    def test09_0(self):
        """NPD Anisotropic Atom [Mitch Miller]"""

        imol_miller = handle_read_draw_molecule_with_recentre(m_miller_res, 1)

        self.failUnless(valid_model_molecule_qm(imol_miller),
                        "Bad read of miller test molecule")

        set_show_aniso(1)  # crash?
        rotate_y_scene(rotate_n_frames(100), 0.1)
        close_molecule(imol_miller)
        set_show_aniso(0)
                    
        
# cheesy test, close the shelx molecules
#             
    def test10_0(self):
        """close shelx molecules"""

        close_molecule(imol_insulin_map)
        close_molecule(imol_insulin_res)

        self.failIf(valid_model_molecule_qm(imol_insulin_res))
        self.failIf(valid_map_molecule_qm(imol_insulin_map))


    def test09_0(self):
        """Aniso Bs in P21"""

        from types import ListType

        def aniso_b_from_atom(atom):
            occ_etc = atom[1]
            return occ_etc[1]

        imol = unittest_pdb("horma-p21.res")
        self.failUnless(valid_model_molecule_qm(imol),
                        "  failed to get a valid imoll from horma-p21.res")
        write_shelx_ins_file(imol, "new-horma.ins")
        imol_2 = read_pdb("new-horma.ins")
        at_1 = get_atom(imol   , "A", 4, "", " N1 ")
        at_2 = get_atom(imol_2 , "A", 4, "", " N1 ")
        b_1  = aniso_b_from_atom(at_1)
        b_2  = aniso_b_from_atom(at_2)

        print "b_1:", b_1
        print "b_2:", b_2

        self.failUnless(type(b_1) is ListType)
        self.failUnless(type(b_2) is ListType)

        map(lambda x, y: self.failUnlessAlmostEqual(x, y), b_1, b_2)
        

        
