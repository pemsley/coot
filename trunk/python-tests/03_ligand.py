import unittest
import os

class LigandTestFunctions(unittest.TestCase):
    
    def test01_0(self):
        """Get monomer test"""

        if (have_test_skip):
            self.skipIf(not have_ccp4_qm, "CCP4 not set up - skipping 3GP test")
        else:
            if (not have_ccp4_qm):
                print "CCP4 not set up - skipping 3GP test (actually passing!)"
                skipped_tests.append("Get monomer test (no CCP4)")
                return
	# BL says:: we shoudl change the monomer_...., so that it can take 2 args
	imol = monomer_molecule_from_3_let_code("3GP", "", "")
	if (valid_model_molecule_qm(imol)):
            global imol_ligand
            imol_ligand = imol 
            delete_residue_hydrogens(imol, "A", 1, "", "")


    def test02_0(self):
	"""Set Bond thickness test"""

        if (have_test_skip):
            self.skipIf(not valid_model_molecule_qm(imol_ligand), "   No ligand molecule - Skipping bond thickness test")
        else:
            if (not valid_model_molecule_qm(imol_ligand)):
                print "   No ligand molecule - Skipping bond thickness test (actually passing!)"
                skipped_tests.append("Set Bond thickness test")
                return
            
	set_bond_thickness(imol_ligand, 5)


    def test03_0(self):
	"""Move and Refine Ligand test"""
        
	new_rc = [55.3, 9.1, 20.6]
        if (have_test_skip):
            self.skipIf(not valid_model_molecule_qm(imol_ligand), "no ligand - skipping test")
        else:
            if (not valid_model_molecule_qm(imol_ligand)):
                print "no ligand (skipping - actually passing test!)"
                skipped_tests.append("Move and Refine Ligand test")
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
        
        # note to self, make sure this file is available when then test
        # is run for real.
        #
        if (have_test_skip):
            self.skipIf(not os.path.isfile("monomer-3GP.pdb"), "pdb file not found - skipping flip residue test")
        else:
            if (not os.path.isfile("monomer-3GP.pdb")):
                print "pdb file does not exist - skipping flip residue test (actually passing!)"
                skipped_tests.append("flip residue (no pdb file)")
                return
        imol_orig = read_pdb("monomer-3GP.pdb")
        if (have_test_skip):
            self.skipIf(not valid_model_molecule_qm(imol_orig), "pdb file not found - skipping flip residue test")
        else:
            if (not valid_model_molecule_qm(imol_orig)):
                print "pdb file not found - skipping flip residue test (actually passing!)"
                skipped_tests.append("flip residue (no pdb file)")
                return
        imol_copy = copy_molecule(imol_orig)
        # BL extra
        # we go to the molecule otherwise we dont pick up the right active_atom
        set_go_to_atom_molecule(imol_orig)
        set_go_to_atom_chain_residue_atom_name("A", 1, " C8 ")
        # end BL etxra
        active_atom = active_residue()
        if (have_test_skip):
            self.skipIf(not active_atom, "No active atom found - skipping flip residue test")
        else:
            if (not active_atom):
                print "No active atom - skipping flip residue test (actually passing!)"
                skipped_tests.append("flip residue test (no active atom)")
                return
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        self.failIf(imol == imol_orig, "oops - didn't pick the copy for active res")
        flip_ligand(imol, chain_id, res_no)
        atom_orig_1 = get_atom(imol_orig, "A", 1, " C8 ")
        atom_move_1 = get_atom(imol     , "A", 1, " C8 ")
        d = bond_length(atom_orig_1[2], atom_move_1[2])
        print "distance: ", d
        self.failUnless(d > 2.1, "fail to move test atom d1")
        flip_ligand(imol, chain_id, res_no)
        flip_ligand(imol, chain_id, res_no)
        flip_ligand(imol, chain_id, res_no)
        #  having flipped it round the axes 4
        # times, we should be back where we
        # started.
        atom_orig_1 = get_atom(imol_orig, "A", 1, " C8 ")
        atom_move_1 = get_atom(imol     , "A", 1, " C8 ")
        d2 = bond_length(atom_orig_1[2], atom_move_1[2])
        print "distance d2: ", d2
        self.failUnless(d2 < 0.001, "fail to move atom back to start d2")


            
