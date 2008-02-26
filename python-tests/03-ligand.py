import unittest

#rnase_pdb = "C:\\MinGW\\msys\\home\\bernhard\\coot-sample\\demo-bad.pdb"
#rnase_mtz = "C:\\MinGW\\msys\\home\\bernhard\\coot-sample\\rnasa-1.8-all_refmac1.mtz"
#terminal-residue-test-pdb = ""

have_ccp4_qm = True
global imol_ligand
imol_ligand = -1


class LigandTestFunctions(unittest.TestCase):
    
    def test11(self):
        """Get Monomer test"""
	self.skipIf(not have_ccp4_qm, "CCP4 not set up - skipping 3GP test")
	# BL says:: we shoudl change the monomer_...., so that it can take 2 args
	imol = monomer_molecule_from_3_let_code("3GP", "", "")
	if (valid_model_molecule_qm(imol)):
		global imol_ligand
		imol_ligand = imol 
		print "BL DEBUG:: imol_lig", imol_ligand
		delete_residue_hydrogens(imol, "A", 1, "", "")

    def test12(self):
	"""Set Bond thickness test"""
	print "BL DEBUG:: imol_lig", imol_ligand
	print "BL DEBUG:: valdi", valid_model_molecule_qm(imol_ligand)
	self.skipIf(not valid_model_molecule_qm(imol_ligand), "No ligand molecule - Skipping bond thickness test")
	set_bond_thickness(imol_ligand, 5)

    def test4(self):
	"""Move and Refine Ligand test"""
	new_rc = [55.3, 9.1, 20.6]
	global imol_ligand
	self.skipIf(not valid_model_molecule_qm(imol_ligand), "no ligand")
	# set the view
	view_number = add_view([54.5698,8.7148,20.5308],
				[0.046229,-0.157139,-0.805581,0.569395],
				19.8858,
				"ligand_view")
# BL says: for now no return on add_view!
#	go_to_view_number(view_number, 1)
	go_to_view_number(1, 1)

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
	rotate_y_scene(600, 0.1)
	if (replacement_state == 0):
		set_refinement_immediate_replacement(0)
	if (backup_mode == 1):
		turn_on_backup(imol_ligand)

