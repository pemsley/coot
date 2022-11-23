# 07_ncs.py
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
import ncs # coot.ncs one day
import begin


def insulin_res():
    return os.path.join(begin.unittest_data_dir, "insulin.res")
def imol_insulin():
    return coot.read_pdb(begin.insulin_res())

class NcsTestFunctions(unittest.TestCase):


    def test01_1(self):
            """NCS maps test"""

            global imol_rnase_map
            global imol_rnase
            self.assertTrue(coot_utils.valid_model_molecule_qm(imol_rnase), "imol_rnase not valid")

            self.assertTrue(coot_utils.valid_map_molecule_qm(imol_rnase_map), "imol_rnase_map not valid")

            n_mols = coot.graphics_n_molecules()
            # try to make it trip up by doing it twice:
            imol_map_2 = coot.make_and_draw_map(begin.rnase_mtz(), "FWT", "PHWT", "", 0 ,0)
            coot.make_dynamically_transformed_ncs_maps(imol_rnase, imol_rnase_map, 0)
            coot.make_dynamically_transformed_ncs_maps(imol_rnase, imol_map_2, 0)
            # 2*2 + 1 new maps should have been made
            n_new = coot.graphics_n_molecules()
            self.assertEqual(n_new, (n_mols + 5), "no match in number of molecules %s %s" %(n_mols, n_new))


    def test02_0(self):
        """NCS chains info"""

        # should return False
        ncs_chain_info = coot.ncs_chain_ids_py(-1)
        self.assertFalse(ncs_chain_info, "   Fail: ncs-chains returns %s, should be False" %ncs_chain_info)

        # a normal case
        coot.make_ncs_ghosts_maybe(imol_rnase)
        ncs_chain_info = coot.ncs_chain_ids_py(imol_rnase)
        self.assertTrue(ncs_chain_info, "   Fail: ncs-chain-ids returns False")
        self.assertTrue(len(ncs_chain_info) > 0, "   Fail: ncs-chains returns %s" %ncs_chain_info)
        first_ghost = ncs_chain_info[0]
        self.assertTrue(len(first_ghost) > 1, "    Fail: first-ghost %s" %first_ghost)
        print("   NCS info: ", ncs_chain_info)


    def test03_0(self):
        """NCS deviation info"""

        # should return False
        ncs_chain_info = coot.ncs_chain_differences(-1, "XX")
        self.assertFalse(ncs_chain_info, "   Fail: ncs-chains returns %s, should be False" %ncs_chain_info)

        # should return False for insulin
        #global imol_insulin
        ncs_chain_info = coot.ncs_chain_differences(imol_insulin(), "A")
        self.assertFalse(ncs_chain_info, "   Fail: ncs-chains for insulin returns %s, should be False" %ncs_chain_info)

        # a normal case
        coot.make_ncs_ghosts_maybe(imol_rnase)
        ncs_chain_info = coot.ncs_chain_differences(imol_rnase, "A")
        self.assertTrue(ncs_chain_info, "   Fail: ncs-chain-differences returns False")
        self.assertEqual(len(ncs_chain_info), 3,
                             """   Fail on length: length ncs-chain-differences should be 3 is %s\n
                                     ncs-chain-differences returns %s""" %(len(ncs_chain_info), ncs_chain_info))


# Phil spotted this bug (it wouldn't refine across the residue 3->4
# peptide bond - because out out of order residues
#
    def test04_0(self):
        """NCS Residue Range copy"""

        # ls is a list of numbers. Is it in ascending order? Return True or False
        def ascending_order_qm(ls):
            asc = sorted(ls)
            if (ls == asc):
                return True
            else:
                return False

        # prepare the input
        imol = coot.read_pdb(begin.rnase_pdb())
        self.assertTrue(coot_utils.valid_model_molecule_qm(imol), "fail to read %s" %begin.rnase_pdb())
        for r in range(1, 4):
            coot.delete_residue(imol, "B", r, "")
        # make ghosts
        coot.make_ncs_ghosts_maybe(imol)
        # evaluate the function to be tested
        coot.copy_residue_range_from_ncs_master_to_others(imol, "A", 1, 6)
        # check the result
        chain_id = "B"
        n_residues = coot.chain_n_residues(chain_id, imol)
        seqnum_order = [coot.seqnum_from_serial_number(imol, chain_id, serial_number) for serial_number in range(n_residues)]

        # return a boolean value
        self.assertTrue(ascending_order_qm(seqnum_order), "fail with seqnum_order %s not ascending" %seqnum_order)

        # Now we can test a mutation.
        #
        # check that the mutate function returns success.
        #
        # check that the new residue type is indeed a TRP
        #
        # make the NCS copy
        #
        # check that the NCS target residue has had it's type
        # updated to TRP.
        #
        mutate_success = coot.mutate(imol, "A", 2, "", "TRP")
        self.assertTrue(mutate_success == 1, "Mutate fails.")

        rname = coot.residue_name(imol, "A", 2, "")
        self.assertTrue(rname == "TRP", "Mutate fails - master not a TRP. %s" %rname)

        coot.copy_residue_range_from_ncs_master_to_others(imol, "A", 1, 3)
        rname = coot.residue_name(imol, "B", 2, "")
        self.assertTrue(rname == "TRP", "Mutate fails - peer not a TRP (%s)" %rname)


# Perhaps combine this with previous test, its tests the same function.
#
    def test05_0(self):
        """NCS Residue Range edit to all chains"""

        imol = begin.unittest_pdb("pdb1t6q.ent")
        self.assertTrue(coot_utils.valid_model_molecule_qm(imol))

        coot.mutate(imol, "A", 50, "", "ASP")
        ncs.skip_to_next_ncs_chain("forward")     # generate the ghosts

        coot.copy_residue_range_from_ncs_master_to_others(imol, "A", 50, 50)

        # did it apply?
        result = []
        for chain_id_peer in ["B", "C"]:
            resname = coot.residue_name(imol, chain_id_peer, 50, "")
            if (resname == "ASP"):
                result.append(True)
            else:
                result.append(False)
        print("result:", result)
        self.assertTrue(all(result))


# This excercises a failure reported by Engin Ozkan 20081209.  Oh
# D'oh, I'd hard-coded "A" as the master chain id into
# manual-ncs-ghosts, I should have used the beginning of the list of
# chain-ids instead.
#
    def test06_0(self):
        """Manual NCS ghosts generates correct NCS chain ids"""

        imol = begin.unittest_pdb("pdb1hvv.ent")

        coot.set_draw_ncs_ghosts(imol, 1)
        coot.ncs_control_change_ncs_master_to_chain_id(imol, "B")
        coot.make_ncs_ghosts_maybe(imol)
        ncs_ghost_chains_1 = coot_utils.ncs_chain_ids(imol)
        coot.manual_ncs_ghosts(imol, 220, 230, ["B", "A", "C", "D"])
        ncs_ghost_chains_2 = coot_utils.ncs_chain_ids(imol)

        print("   NCS ghost chain IDs pre:  ", ncs_ghost_chains_1)
        print("   NCS ghost chain IDs post: ", ncs_ghost_chains_2)

        self.assertEqual(ncs_ghost_chains_1, [["B", "A", "C", "D"]])
        self.assertEqual(ncs_ghost_chains_2, [["B", "A", "C", "D"]])


    def test07_0(self):
        """NCS maps overwrite existing maps"""

        # first close all the maps that have "NCS found" in the name:
        for imol in coot_utils.molecule_number_list():
            if ("NCS found" in coot.molecule_name(imol)):
                coot.close_molecule(imol)

        imol = begin.unittest_pdb("pdb1hvv.ent")
        imol_map = coot.make_and_draw_map(os.path.join(begin.unittest_data_dir, "1hvv_sigmaa.mtz"),
                                          "2FOFCWT", "PH2FOFCWT", "", 0, 0)

        coot.make_dynamically_transformed_ncs_maps(imol, imol_map, 0)
        coot.make_dynamically_transformed_ncs_maps(imol, imol_map, 0)
        coot.make_dynamically_transformed_ncs_maps(imol, imol_map, 1)

        result_list = []
        molecule_names = list(map(coot.molecule_name, coot_utils.molecule_number_list()))
        print("BL DEBUG:: molecule_names", molecule_names)
        for chain_id in ["B", "C", "D"]:
            test_name = "Map " + \
                        str(imol_map) + " " + \
                        "NCS found from matching Chain " + \
                        chain_id + \
                        " onto Chain A"

            n_matchers = molecule_names.count(test_name)
            print("BL DEBUG:: n_matchers", n_matchers)
            self.assertTrue(n_matchers >= 2, "  Failed to find matching NCS chain %s" %chain_id)
