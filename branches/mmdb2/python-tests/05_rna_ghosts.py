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
            for chain_id in chain_ids(imol):
                if (not is_solvent_chain_qm(imol, chain_id)):
                    n_residues = chain_n_residues(chain_id, imol)
                    print "   There are %s residues in chain %s" %(n_residues, chain_id)

                    for serial_number in range(n_residues):

                        res_name = resname_from_serial_number(imol, chain_id, serial_number)
                        res_no   = seqnum_from_serial_number (imol, chain_id, serial_number)
                        ins_code = insertion_code_from_serial_number(imol, chain_id, serial_number)
                        atom_ls  = residue_info(imol, chain_id, res_no, ins_code)
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
            set_atom_attributes(atom_attribute_settings)

        # main body
        rna_mol = ideal_nucleic_acid("RNA", "A", 0, "GACUCUAG")
        copy_rna_mol = copy_molecule(rna_mol)

        # move the view over a bit so we can see the atoms being jiggled
        rc = rotation_centre()
        set_rotation_centre(rc[0] + 12,
                            rc[1] + 3,
                            rc[2])
        view_number = add_view([74.7079, 10.6267, 24.3308],
                               [-0.713385, -0.0433099, -0.105865, -0.691373],
                               70.3919,
                               "RNA-builder-view")
        go_to_view_number(view_number, 1)

        # now jiggle the atoms of copy-rna-mol
        jiggle_atoms_of_mol(copy_rna_mol)
        merge_molecules([copy_rna_mol], rna_mol)

        imol_copy = copy_molecule(rna_mol)
        clear_lsq_matches()
        add_lsq_match(1, 6, "A", 1, 6, "C", 0)   # ref mov - all atoms
        rrtop = apply_lsq_matches(imol_copy, imol_copy)
        rtop = []
        for r in rrtop:
            rtop.extend(r)
        #close_molecule(imol_copy)
        #close_molecule(copy_rna_mol)
        self.failUnless(rtop, "Failed to get matching matrix")
        set_draw_ncs_ghosts(rna_mol, 1)
        add_ncs_matrix(rna_mol, "C", "A", *rtop)
        view_number = add_view([72.3306, 10.6899, 24.073],
                               [-0.240736, -0.674651, -0.690658, -0.0994136],
                               14.9021,
                               "RNA-ghots-view")
        go_to_view_number(view_number, 1)
        rotate_y_scene(rotate_n_frames(200), 1)
