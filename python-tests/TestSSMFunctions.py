# 06_ssm.py
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

ssm_ref_1_file = os.path.join(coot_testing_utils.unittest_data_dir, "2qd9.pdb")
ssm_ref_2_file = os.path.join(coot_testing_utils.unittest_data_dir, "2gtn.pdb")

class TestSSMFunctions(unittest.TestCase):

    def test01_0(self):
        """SSM - Frank von Delft's Example"""

        imol_a = coot.handle_read_draw_molecule_with_recentre(os.path.join(coot_testing_utils.unittest_data_dir, "1wly.pdb"), 0)
        imol_b = coot.handle_read_draw_molecule_with_recentre(os.path.join(coot_testing_utils.unittest_data_dir, "1yb5.pdb"), 1)

        self.assertTrue(coot_utils.valid_model_molecule_qm(imol_a) and coot_utils.valid_model_molecule_qm(imol_b))

        coot.graphics_to_ca_plus_ligands_representation(imol_a)

        coot.superpose_with_atom_selection(imol_a, imol_b, "A/2-111", "A/6-115", 0)
        coot.set_rotation_centre(65.65, -3, -4)
        view_number = coot_utils.add_view([49.7269, 7.69693, 3.93221],
                                          [-0.772277, 0.277494, 0.292497, 0.490948],
                                          98.9608, "SSM View")
        coot.go_to_view_number(view_number, 1)
        # coot.rotate_y_scene(coot.rotate_n_frames(100), 0.1)
        coot.set_mol_displayed(imol_a, 0)
        coot.set_mol_displayed(imol_b, 0)
        # didnt crash....

    def test02_0(self):
        """SSM - Alice Dawson's Example"""
        imol_s = coot.handle_read_draw_molecule_with_recentre(os.path.join(coot_testing_utils.unittest_data_dir, "1pyd.pdb"), 0)

        coot.graphics_to_ca_plus_ligands_representation(imol_s)
        coot.set_graphics_window_size(678, 452)

        coot_utils.print_molecule_names()

        coot.superpose_with_atom_selection(imol_s, imol_s, "A/100-400", "B/50-450", 1)
        imol_copy = coot.graphics_n_molecules() - 1
        coot.graphics_to_ca_plus_ligands_representation(imol_copy)
        # coot.rotate_y_scene(coot.rotate_n_frames(100), 0.1)
        # didnt crash...

    def test03_0(self):
        """SSM by atom selection [JED Example]"""

        imol_1 = coot.read_pdb(ssm_ref_1_file)
        imol_2 = coot.read_pdb(ssm_ref_2_file)
        coot.superpose_with_atom_selection(imol_1, imol_2,"//A/140-160", "//A/140-160", 0)
        # didn't crash?
        # testing for something else?!

