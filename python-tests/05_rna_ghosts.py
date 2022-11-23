# 05_rna_ghosts.py
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

class RnaGhostsTestFunctions(unittest.TestCase):

    def test01_0(self):
        """RNA NCS Ghosts"""
        import random

        jiggle_random = random.uniform(-0.5, 0.5)

        # a bit of ugliness here, we set! additions to
        # atom-attribute-settings, rather than consing them with
        # recursion over chain-id, serial-number and atom-list.
        # Hmm...  Come back to improve one day?
        #
        def jiggle_atoms_of_mol(imol):
            atom_attribute_settings = []
            for chain_id in coot_utils.chain_ids(imol):
                if (not coot_utils.is_solvent_chain_qm(imol, chain_id)):
                    n_residues = coot.chain_n_residues(chain_id, imol)
                    print("   There are %s residues in chain %s" %(n_residues, chain_id))

                    for serial_number in range(n_residues):

                        res_name = coot.resname_from_serial_number(imol, chain_id, serial_number)
                        res_no   = coot.seqnum_from_serial_number (imol, chain_id, serial_number)
                        ins_code = coot.insertion_code_from_serial_number(imol, chain_id, serial_number)
                        atom_ls  = coot.residue_info(imol, chain_id, res_no, ins_code)
                        for atom in atom_ls:
                            compound_name = atom[0]
                            atom_name = compound_name[0]
                            alt_conf  = compound_name[1]
                            xyz = atom[2]
                            x = xyz[0]
                            y = xyz[1]
                            z = xyz[2]
                            atom_attribute_settings.append([imol, chain_id, res_no, ins_code, atom_name, alt_conf,
                                                             "x", (24 + x + 0.3 * jiggle_random)])
                            atom_attribute_settings.append([imol, chain_id, res_no, ins_code, atom_name, alt_conf,
                                                             "y", (24 + y + 0.3 * jiggle_random)])
                            atom_attribute_settings.append([imol, chain_id, res_no, ins_code, atom_name, alt_conf,
                                                             "z", (24 + z + 0.3 * jiggle_random)])
                            coot.set_atom_attributes(atom_attribute_settings)

        # main body
        rna_mol = coot.ideal_nucleic_acid("RNA", "A", 0, "GACUCUAG")
        copy_rna_mol = coot.copy_molecule(rna_mol)

        # move the view over a bit so we can see the atoms being jiggled
        rc = coot.rotation_centre()
        coot.set_rotation_centre(rc[0] + 12,
                                 rc[1] + 3,
                                 rc[2])
        view_number = coot.add_view([74.7079, 10.6267, 24.3308],
                                    [-0.713385, -0.0433099, -0.105865, -0.691373],
                                    70.3919,
                                    "RNA-builder-view")
        coot.go_to_view_number(view_number, 1)

        # now jiggle the atoms of copy-rna-mol
        coot.jiggle_atoms_of_mol(copy_rna_mol)
        coot.merge_molecules([copy_rna_mol], rna_mol)

        imol_copy = coot.copy_molecule(rna_mol)
        coot.clear_lsq_matches()
        coot.add_lsq_match(1, 6, "A", 1, 6, "C", 0)   # ref mov - all atoms
        rrtop = coot.apply_lsq_matches(imol_copy, imol_copy)
        rtop = []
        for r in rrtop:
            rtop.extend(r)
        #close_molecule(imol_copy)
        #close_molecule(copy_rna_mol)
        self.assertTrue(rtop, "Failed to get matching matrix")
        coot.set_draw_ncs_ghosts(rna_mol, 1)
        coot.add_ncs_matrix(rna_mol, "C", "A", *rtop)
        view_number = coot.add_view([72.3306, 10.6899, 24.073],
                                    [-0.240736, -0.674651, -0.690658, -0.0994136],
                                    14.9021,
                                    "RNA-ghots-view")
        coot.go_to_view_number(view_number, 1)
        coot.rotate_y_scene(coot.rotate_n_frames(200), 1)
