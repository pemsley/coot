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

imol_hof_res = False

insulin_fcf = os.path.join(unittest_data_dir, "insulin.fcf")
insulin_res = os.path.join(unittest_data_dir, "insulin.res")
hollander_ins = os.path.join(unittest_data_dir, "hollander.ins")
global imol_insulin_res
imol_insulin_res = -1 # see later
m_miller_res = os.path.join(unittest_data_dir, "miller/shelx-test4-NPD-mini.res")


class ShelxTestFunctions(unittest.TestCase):

    def test01_0(self):
        """Read small molecule .res file"""

        if (have_test_skip):
            self.skipIf(not type(hof_res) is StringType, "hof-res not defined - skipping test")
            self.skipIf(not os.path.isfile(hof_res), "%s does not exist - skipping test" %hof_res)
        else:
            if (not type(hof_res) is StringType):
                print "hof-res not defined - skipping test (actually passing!)"
                skipped_tests.append("Read small molecule .res file")
                return
            if (not os.path.isfile(hof_res)):
                print "%s does not exist - skipping test (actually passing!)" %hof_res
                skipped_tests.append("Read small molecule .res file")
                return
            
        imol = read_pdb(hof_res)
        self.failUnless(valid_model_molecule_qm(imol))
        imol_hof_res = imol


    def test02_0(self):
        """Read hollander small molecule .res file"""

        if (have_test_skip):
            self.skipIf(not os.path.isfile(hollander_ins), "%s  does not exist - skipping test" %hollander_ins)
        else:
            if (not os.path.isfile(hollander_ins)):
                print "%s  does not exist - skipping test (actually passing)" %hollander_ins
                skipped_tests.append("Read hollander small molecule .res file")
                return

        imol = read_pdb(hollander_ins)
        self.failUnless(valid_model_molecule_qm(imol), "   fail: bad molecule for %s" %hollander_ins)
        spg = show_spacegroup(imol)
        self.failUnlessEqual(spg, "I 41 2 2", "   fail: wrong spacegroup for %s %s" %(hollander_ins, spg))
        

    def test03_0(self):
        """read shelx insulin with fcf"""

        global imol_insulin_res
        
        imol_insulin_res_local = handle_read_draw_molecule_with_recentre(insulin_res, 1)
        self.failUnless(valid_model_molecule_qm(imol_insulin_res_local), "   Bad insulin.res: %s for %s" %(insulin_res, imol_insulin_res_local))

        imol_insulin_res = imol_insulin_res_local  # used in water addition test
        imol = handle_shelx_fcf_file(insulin_fcf)
        self.failUnless(valid_map_molecule_qm(imol), "    Bad read of %s %s" %(insulin_fcf, imol))
        name = molecule_name(imol)
        cif_name = insulin_fcf + ".cif SigmaA"
        self.failUnlessEqual(name, cif_name, "   Bad name match %s != %s" %(name, cif_name))

        self.failUnlessEqual(show_spacegroup(imol), show_spacegroup(imol_insulin_res_local),
                             "   Mismatch spacegroups %s %s" %(show_spacegroup(imol),
                                                               show_spacegroup(imol_insulin_res_local)))

        self.failUnlessEqual(show_spacegroup(imol), "I 21 3", "   Bad spacegroups %s" %(show_spacegroup(imol)))

        # good then
        rotate_y_scene(rotate_n_frames(200), 0.1)
        close_molecule(imol)
        # imol_insulin_res_local needed later -> global?
        

    # The "SHELX Woe" problem
    #
    def test04_0(self):
        """Write an INS from PDB test"""

        # First, check return status on a bogus molecule
        status = write_shelx_ins_file(204050, "no_molecule.ins")
        self.failUnlessEqual(status, 0, "bad exit status from write_shelx_ins_file on bogus molecule %s" %status)

        # Now check return status on a bogus molecule.  The Happy Path.
        #
        self.failUnless(valid_model_molecule_qm(imol_rnase), "   imol-rnase not valid.")
        rnase_ins = "rnase.ins"
        status = write_shelx_ins_file(imol_rnase, rnase_ins)
        self.failUnlessEqual(status, 1, "   failure to write INS file %s from PDB: status %s" %(rnase_ins, status))


    def test05_0(self):
        """Add water to SHELX molecule"""

        global imol_insulin_res

        set_pointer_atom_molecule(imol_insulin_res)
        set_rotation_centre(3, -1, 60)
        place_typed_atom_at_pointer("Water")
        # test is to have to occupancy of the new HOH to be 11.0
        chain_id = water_chain(imol_insulin_res)
        n_residues = chain_n_residues(chain_id, imol_insulin_res)
        serial_number = n_residues - 1
        res_name = resname_from_serial_number(imol_insulin_res, chain_id, serial_number)
        res_no   = seqnum_from_serial_number (imol_insulin_res, chain_id, serial_number)
        ins_code = insertion_code_from_serial_number(imol_insulin_res, chain_id, serial_number)

        self.failUnlessEqual(res_name, "HOH")
        atom_list = residue_info(imol_insulin_res, chain_id, res_no, ins_code)
        self.failIfEqual(atom_list, [])
        for atom in atom_list:
            occ = atom[1][0]
            self.failUnlessAlmostEqual(occ, 11.0, 1, "  bad occupancy in SHELXL molecule %s" %atom)


# non positive definite anistropic atom (reported by Mitch Miller)
# crash test
    def test06_0(self):
        """NPD Anisotripic Atom [Mitch Miller]"""

        imol_miller = handle_read_draw_molecule_with_recentre(m_miller_res, 1)

        self.failUnless(valid_model_molecule_qm(imol_miller),
                        "Bad read of miller test molecule")

        set_show_aniso(1)  # crash?
        rotate_y_scene(rotate_n_frames(100), 0.1)
        close_molecule(imol_miller)
        set_show_aniso(0)
                    
        
            
