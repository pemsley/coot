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

        for i in range(0,50):
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
        
        
