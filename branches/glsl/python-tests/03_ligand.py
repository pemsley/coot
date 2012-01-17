#    03_ligand.oy
#    Copyright (C) 2008  Bernhard Lohkamp, The University of York
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import unittest
import os

class LigandTestFunctions(unittest.TestCase):
    
    def test01_0(self):
        """Get monomer test"""

        if (not have_ccp4_qm):
            print "No CCP4 - Copying in test files"
            for file_name in ["monomer-3GP.pdb", "libcheck_3GP.cif"]:
                f_full = os.path.join(unittest_data_dir, file_name)
                t_file = os.path.join("coot-ccp4", file_name)
                if (os.path.isfile(f_full)):
                    import shutil
                    make_directory_maybe("coot-ccp4")
                    print "   copy-file", f_full, t_file
                    shutil.copyfile(f_full, t_file)
                else:
                    if self.skip_test(True, "Cannot find file %s, skipping test" %f_full):
                        return
                
	imol = monomer_molecule_from_3_let_code("3GP", "")
	if (valid_model_molecule_qm(imol)):
            global imol_ligand
            imol_ligand = imol 
            delete_residue_hydrogens(imol, "A", 1, "", "")


    def test02_0(self):
	"""Set Bond thickness test"""

        if self.skip_test(not valid_model_molecule_qm(imol_ligand),
                          "   No ligand molecule - Skipping bond thickness test"):
            return
            
	set_bond_thickness(imol_ligand, 5)


    def test02_1(self):
        """Delete all-molecule Hydrogens"""
        imol = unittest_pdb("monomer-3GP.pdb")
        self.failUnless(valid_model_molecule_qm(imol))
        print "here 1"
        n = delete_hydrogens(imol)
        self.failUnless(n > 0)
        

    def test03_0(self):
	"""Move and Refine Ligand test"""
        
	new_rc = [55.3, 9.1, 20.6]
        if self.skip_test(not valid_model_molecule_qm(imol_ligand),
                          "no ligand - skipping test"):
            return

	# set the view
	view_number = add_view([54.5698,8.7148,20.5308],
				[0.046229,-0.157139,-0.805581,0.569395],
				19.8858,
				"ligand_view")
	go_to_view_number(view_number, 1)

	# update the map:
	set_rotation_centre(*new_rc)
	move_molecule_here(imol_ligand)
	backup_mode = backup_state(imol_ligand)
	alt_conf = ""
	replacement_state = refinement_immediate_replacement_state()
        
	turn_off_backup(imol_ligand)
	set_refinement_immediate_replacement(1)
	refine_zone(imol_ligand, "A", 1, 1, alt_conf)
	accept_regularizement()
	rotate_y_scene(rotate_n_frames(600), 0.1)
	if (replacement_state == 0):
		set_refinement_immediate_replacement(0)
	if (backup_mode == 1):
		turn_on_backup(imol_ligand)

        # testing what?

    def test04_0(self):
        """Many Molecules - Ligand Fitting"""

        npo_pdb = os.path.join(unittest_data_dir, "monomer-NPO.pdb")
        pdb43ca_pdb = os.path.join(unittest_data_dir, "pdb43ca-sans-NPO-refmaced.pdb")
        pdb43ca_mtz = os.path.join(unittest_data_dir, "pdb43ca-sans-NPO-refmaced.mtz")
        imol_npo = handle_read_draw_molecule_with_recentre(npo_pdb, 0)

        self.failUnless(valid_model_molecule_qm(imol_npo), "no valid molecule")

        for i in range(0, 5):
            imol_copy = copy_molecule(imol_npo)
            set_mol_displayed(imol_copy, 0)

        imol_protein = read_pdb(pdb43ca_pdb)
        imol_map_2 = auto_read_make_and_draw_maps(pdb43ca_mtz)
        imol_map_1 = imol_map_2 - 1

        add_ligand_clear_ligands()
        set_ligand_search_protein_molecule(imol_protein)
        set_ligand_search_map_molecule(imol_map_1)
        add_ligand_search_ligand_molecule(imol_npo)

        solutions = execute_ligand_search()   # crash?
        print "   Fitting NPO gave these results: ", solutions
        set_map_displayed(imol_map_1, 0)
        set_map_displayed(imol_map_2, 0)
        set_mol_displayed(imol_protein, 0)
        set_mol_displayed(imol_npo, 0)
        # checked for non crash

        
    def test05_0(self):
        """flip residue (around eigen vectors)"""

        mon_file = os.path.join("coot-ccp4", "monomer-3GP.pdb")
        self.failIf(not os.path.isfile(mon_file),
                    "  Oops! file not found! coot-ccp4/monomer-3GP.pdb")

        imol_orig = read_pdb(mon_file)
        imol_copy = copy_molecule(imol_orig)

        self.failIf(not valid_model_molecule_qm(imol_orig), "not valid molecule for monomer-3GP.pdb")

        self.failIf(not valid_model_molecule_qm(imol_copy), "not valid molecule for copy of monomer-3GP.pdb")

        # BL extra
        # we go to the molecule otherwise we dont pick up the right active_atom
        set_go_to_atom_molecule(imol_orig)
        set_go_to_atom_chain_residue_atom_name("A", 1, " C8 ")
        # end BL etxra
        active_atom = active_residue()
        if self.skip_test(not active_atom,
                          "No active atom found - skipping flip residue test"):
            return
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        self.failIf(imol == imol_orig, "oops - didn't pick the copy for active res")
        flip_ligand(imol, chain_id, res_no)
        atom_orig_1 = get_atom(imol_orig, "A", 1, "", " C8 ")
        atom_move_1 = get_atom(imol     , "A", 1, "", " C8 ")

        from types import ListType

        self.failUnless(type(atom_orig_1) is ListType, "atom_orig_1 not found")

        self.failUnless(type(atom_move_1) is ListType, "atom_move_1 not found")

        d = bond_length(atom_orig_1[2], atom_move_1[2])
        print "distance: ", d
        self.failUnless(d > 2.1, "fail to move test atom d1")
        flip_ligand(imol, chain_id, res_no)
        flip_ligand(imol, chain_id, res_no)
        flip_ligand(imol, chain_id, res_no)
        #  having flipped it round the axes 4
        # times, we should be back where we
        # started.
        atom_orig_1 = get_atom(imol_orig, "A", 1, "", " C8 ")
        atom_move_1 = get_atom(imol     , "A", 1, "", " C8 ")
        d2 = bond_length(atom_orig_1[2], atom_move_1[2])
        print "distance d2: ", d2
        self.failUnless(d2 < 0.001, "fail to move atom back to start d2")


    def test06_0(self):
        """Test dipole"""

        imol = unittest_pdb("dipole-residues.pdb")

        self.failUnless(valid_model_molecule_qm(imol), "dipole-residues.pdb not found")

        residue_specs = [["A", 1, ""],
                         ["A", 2, ""],
                         ["A", 3, ""]]
        dipole = add_dipole_for_residues(imol, residue_specs)

        self.failIf(not dipole, "bad dipole %s" %dipole)

        d = dipole[0]
        dip = dipole[1]

        dip_x = dip[0]
        dip_y = dip[1]
        dip_z = dip[2]
        
        print "info:: dipole components", dip

        self.failUnlessAlmostEqual(dip_y, 0.0, 2, "bad dipole y component %s" %dip_y) 
        self.failUnlessAlmostEqual(dip_z, 0.0, 2, "bad dipole z component %s" %dip_z)

        self.failUnless(dip_x < 0 and dip_x > -20)

        
    def test06_0(self):
        """Reading new dictionary restraints replaces"""

        def get_torsions(r):
            return r["_chem_comp_tor"]

        read_cif_dictionary(os.path.join(unittest_data_dir, "libcheck_3GP.cif"))
        read_cif_dictionary(os.path.join(unittest_data_dir, "libcheck_3GP.cif"))
        read_cif_dictionary(os.path.join(unittest_data_dir, "libcheck_3GP.cif"))

        r = monomer_restraints("3GP")
        self.failUnless(r, "Failed to get restraints from monomer 3GP")
        t = get_torsions(r)

        #
        self.failUnless(len(t) < 26, "torsions: %s %s" %(len(t), t))
        # 22 in new dictionary, it seems
