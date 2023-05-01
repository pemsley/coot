# 01_pdb_mtz.py
# Copyright 2007, 2008 by The University of York
# Copyright 2009, 2010, 2011 by Bernhard Lohkamp
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
# and test that gobject is in place.  
#import gobject

global terminal_residue_test_pdb
terminal_residue_test_pdb = os.path.join(unittest_data_dir,
                                         "tutorial-add-terminal-1-test.pdb")
base_imol = graphics_n_molecules()
rnase_seq = os.path.join(unittest_data_dir, "rnase.seq")

global have_ccp4_qm
have_ccp4_qm = False
global imol_rnase
imol_rnase = -1
global imol_rnase_map
imol_rnase_map = -1
global imol_ligand
imol_ligand = -1
global imol_terminal_residue_test
imol_terminal_residue_test = -1

horne_works_cif   = os.path.join(unittest_data_dir, "lib-B3A.cif")
horne_cif         = os.path.join(unittest_data_dir, "lib-both.cif")
horne_pdb         = os.path.join(unittest_data_dir, "coords-B3A.pdb")
global ins_code_frag_pdb
ins_code_frag_pdb = os.path.join(unittest_data_dir, "ins-code-fragment-pre.pdb")

set_map_radius(4.5) # faster

# CCP4 is set up? If so, set have-ccp4? True
try:
    #global have_ccp4_qm
    ccp4_master = os.getenv("CCP4_MASTER")
    if (os.path.isdir(ccp4_master)):
        print "==== It seems CCP4_MASTER is setup ==="
        print "==== CCP4_MASTER:", ccp4_master
        have_ccp4_qm = True
except:
    print "BL INFO:: Dont have CCP4 master"

class PdbMtzTestFunctions(unittest.TestCase):

    # tests are executed alphanumerical, so we shall give them number,
    # rather than names. We add a 0 in the end to give space for later
    # addition
    def test00_0(self):
        """Post Go To Atom no molecule"""
        post_go_to_atom_window()


    def test01_0(self):
        """Close bad molecule"""
        close_molecule(-2)


    def test02_0(self):
        """Read coordinates test"""
        global imol_rnase
        imol = read_pdb(rnase_pdb())
        imol_rnase = imol
        self.failUnless(valid_model_molecule_qm(imol))


    def test03_0(self):
        """New molecule from bogus molecule"""
        pre_n_molecules = graphics_n_molecules()
        new_mol = new_molecule_by_atom_selection(-5, "//A/0")
        post_n_molecules = graphics_n_molecules()
        self.failUnlessEqual(new_mol, -1 ,
                             "fail on non-matching n molecules (-1)")
        self.failUnlessEqual(pre_n_molecules, post_n_molecules,
                             "fail on non-matching n molecules (=)")


    def test03_1(self):
        """Don't crash on empty NCS from mmCIF file"""

        imol = unittest_pdb("2WF6.pdb")
        print "   closing molecule number", imol
        close_molecule(imol)
         # no testing, just not crashing
         
        
    def test04_0(self):
        """New molecule from bogus atom selection"""

        global imol_rnase
        pre_n_molecules = graphics_n_molecules()
        new_molecule = new_molecule_by_atom_selection(imol_rnase, "//A/100")
        post_n_molecules = graphics_n_molecules()
        # what should happen if there are no atoms in the new-mol?
        print "   INFO:: pre_n_molecules %s   post_n_molecules %s" %(pre_n_molecules, post_n_molecules)
        self.failUnlessEqual(new_molecule, -1,
                             "fail on non-matching n molecules (-1)")
        self.failUnlessEqual(pre_n_molecules, post_n_molecules,
                             "fail on non-matching n molecules (=)")


    def test05_0(self):
        """ins code change and Goto atom over an ins code break"""

        def matches_attributes(atts_obs, atts_ref):
            return atts_obs == atts_ref

        # main line
        self.failUnless(os.path.isfile(ins_code_frag_pdb),
                        "   WARNING:: file not found: %s" %ins_code_frag_pdb)
        frag_pdb = handle_read_draw_molecule_with_recentre(ins_code_frag_pdb, 0)
        set_go_to_atom_molecule(frag_pdb)
        set_go_to_atom_chain_residue_atom_name("A", 68, " CA ")
        ar_1 = active_residue()
        ins_1 = ar_1[3]
        change_residue_number(frag_pdb, "A", 68, "", 68, "A")
        change_residue_number(frag_pdb, "A", 69, "", 68, "B")
        change_residue_number(frag_pdb, "A", 67, "", 68, "")
        ar_2 = active_residue()
        ins_2 = ar_2[3]
        print "   pre and post ins codes: ", ins_1, ins_2
        self.failUnlessEqual(ins_2, "A",
                             "Fail ins code set: %s is not 'A'" %ins_2)

        write_pdb_file(frag_pdb, "post-ins-change.pdb")

        # note note: 68B -> 70 doesn't happen unless we are on 68B (rc distance check)
        #
        test_expected_results = [[goto_next_atom_maybe("A", 67, "",  " CA "),
                                  ["A", 68, "",  " CA "]],
                                 [goto_next_atom_maybe("A", 68, "A", " CA "),
                                  ["A", 68, "B", " CA "]],
                                 [goto_next_atom_maybe("A", 68, "B", " CA "),
                                  ["A", 68, "B",  " CA "]],
                                 [goto_prev_atom_maybe("A", 70, "",  " CA "),
                                  ["A", 68, "B", " CA "]],
                                 [goto_prev_atom_maybe("A", 68, "B", " CA "),
                                  ["A", 68, "A", " CA "]],
                                 [goto_prev_atom_maybe("A", 68, "A", " CA "),
                                  ["A", 68, "",  " CA "]],
                                 [goto_prev_atom_maybe("A", 68, "",  " CA "),
                                  ["A", 66, "",  " CA "]]]

        # Note:: these are taken out beause goto-next-atom-maybe now checks
        # where the screen centre is before the move, so the previously
        # expected results no longer apply.  We could recentre and in these
        # tests - (todo...).
        #[goto_next_atom_maybe("A", 68, "B", " CA "),
        # ["A", 70, "",  " CA "]],
        #[goto_next_atom_maybe("D", 10, "",  " O  "),
        # ["A", 62, "",  " CA "]],

        for test_item in test_expected_results:
            real_result     = test_item[0]
            expected_result = test_item[1]

            self.failUnless(real_result)
            self.failUnlessEqual(real_result,
                                 expected_result,
                                 "   fail: real: %s expected: %s" %(real_result, expected_result))
            print "   pass: ", expected_result


    def test05_1(self):
        """Replace Residue gets correct residue number"""

        imol = unittest_pdb("tutorial-modern.pdb")

        mutate_by_overlap(imol, "A", 86, "PTR")
        rn = residue_name(imol, "A", 86, "")
        self.failUnless(isinstance(rn, str))
        self.failUnless(rn == "PTR")
        # OK, did the the refinement run OK? Check the C-N distance
        N_atom = get_atom(imol, "A", 86, "", " N  ", "")
        C_atom = get_atom(imol, "A", 85, "", " C  ", "")

        dd = bond_length_from_atoms(N_atom, C_atom)
        self.failIf(dd > 1.4)
        self.failIf(dd < 1.25)
        print "C-N dist good enough:", dd


    def test06_0(self):
        """Read a bogus map"""
        pre_n_molecules = graphics_n_molecules()
        imol = handle_read_ccp4_map("bogus.map", 0)
        self.failUnlessEqual(imol, -1,
                             "bogus ccp4 map returns wrong molecule number")
        now_n_molecules = graphics_n_molecules()
        self.failUnlessEqual(now_n_molecules, pre_n_molecules,
                             "bogus ccp4 map creates extra map %s %s " %(pre_n_molecules, now_n_molecules))


    def test07_0(self):
        """Read MTZ test"""

        global imol_rnase_map
        # bogus map test
        pre_n_molecules = graphics_n_molecules()
        imol_map = make_and_draw_map("bogus.mtz", "FWT", "PHWT", "", 5, 6)
        self.failUnlessEqual(imol_map, -1,
                             "   bogus MTZ returns wrong molecule number")
        now_n_molecules = graphics_n_molecules()
        self.failUnlessEqual(now_n_molecules, pre_n_molecules,
                             "   bogus MTZ creates extra map %s %s" %(pre_n_molecules, now_n_molecules))

        # correct mtz test
        imol_map = make_and_draw_map(rnase_mtz(), "FWT","PHWT","",0,0)
        change_contour_level(0)
        change_contour_level(0)
        change_contour_level(0)
        set_imol_refinement_map(imol_map)
        imol_rnase_map = imol_map
        self.failUnless(valid_map_molecule_qm(imol_map))


    def test07_1(self):
        """Auto-read bad MTZ test"""

        mtz_list = ["xx-missing.mtz",
                    os.path.join(unittest_data_dir, "broken.mtz")]

        for file_name in mtz_list:
            r = auto_read_make_and_draw_maps(file_name)
            print "   got status: ", r
        # no crash
        
        
    def test08_0(self):
        """Map Sigma """

        global imol_rnase_map
        global imol_rnase
        self.failUnless(valid_map_molecule_qm(imol_rnase_map))
        v = map_sigma(imol_rnase_map)
        self.failUnless(v > 0.2 and v < 1.0)
        v2 = map_sigma(imol_rnase)
        print "   INFO:: map sigmas", v, v2
        self.failIf(v2)


    def test09_0(self):
        """Another Level Test"""
        imol_map_2 = another_level()
        self.failUnless(valid_map_molecule_qm(imol_map_2))


    def test09_1(self):
        """Sharpen map from map"""
        # don't crash
        mtz_file_name = os.path.join(unittest_data_dir, "3hfl_sigmaa.mtz")
        imol_map = make_and_draw_map(mtz_file_name,
                                     "2FOFCWT", "PH2FOFCWT", "", 0, 0)

        self.failUnless(valid_map_molecule_qm(imol_map),
                        "fail to get map from 3hfl_sigmaa.mtz")

        export_map(imol_map, "test-3hfl.map")
        new_mol = handle_read_ccp4_map("test-3hfl.map", 0)
        self.failUnless(valid_map_molecule_qm(new_mol),
                        "fail to get map from 3hfl map")
        sharpen(new_mol, 5.0)
        # didn't crash


    def test09_1(self):
        """db-main makes mainchain"""
        read_pdb(".")
        imol = read_pdb(rnase_pdb())
        db_mainchain(imol, "A", 10, 20, "forward")
        # didn't hang


    def test09_2(self):
        """"Negative Residues in db-mainchain don't cause a crash"""

        # Oliver Clarke spotted this bug

        imol = unittest_pdb("tutorial-modern.pdb")
        self.failUnless(valid_model_molecule_qm(imol))
        renumber_residue_range(imol, "A", 1, 50, -20)
        imol_mc_1 = db_mainchain(imol, "A", -19, 6, "forwards")
        imol_mc_2 = db_mainchain(imol, "B", 10, 30, "backwards")
        # didnt crash...

    def test10_0(self):
        """Set Atom Atribute Test"""
        atom_ls = []
        global imol_rnase
        set_atom_attribute(imol_rnase, "A", 11, "", " CA ", "", "x", 64.5) # an Angstrom or so
        atom_ls = residue_info(imol_rnase, "A", 11, "")
        self.failIfEqual(atom_ls, [])
        atom = atom_ls[0]
        compound_name = atom[0]
        atom_name = compound_name[0]
        if (atom_name == " CA "):
            x = atom_ls[0][2][0]
            self.failUnlessAlmostEqual(x, 64.5)


    def test11_0(self):
	"""Add Terminal Residue Test"""
	import types

	if (recentre_on_read_pdb() == 0):
            set_recentre_on_read_pdb(1)
        if self.skip_test(not type(terminal_residue_test_pdb) is StringType,
                          "%s does not exist - skipping test" %terminal_residue_test_pdb):
            return
        
        # OK, file exists
        imol = read_pdb(terminal_residue_test_pdb)
        if self.skip_test(not valid_model_molecule_qm(imol), 
                          "%s bad pdb read - skipping test" %terminal_residue_test_pdb):
            return

	set_default_temperature_factor_for_new_atoms(45)
	add_terminal_residue(imol, "A", 1, "ALA", 1)
	write_pdb_file(imol, "regression-test-terminal-residue.pdb")

	# lets modify the test, so that we additionally can test if the added
	# residue is actually in the pdb file we just wrote.
	imol_added = read_pdb("regression-test-terminal-residue.pdb")

	# where did that extra residue go?
	# 
	# Let's copy a fragment from imol and go to
	# the centre of that fragment, and check where
	# we are and compare it to where we expected
	# to be.
	#
	# we shall make sure that we recentre on read pdb (may be overwritten
	# in preferences!!!
	new_mol = new_molecule_by_atom_selection(imol_added, "//A/0")
	self.failUnless(valid_model_molecule_qm(new_mol),
                        "Added residue is not found in the new pdb file")
	move_molecule_here(new_mol)
	rc = rotation_centre()
	ls = [45.6, 15.8, 11.8]
	r = sum([rc[i] - ls[i] for i in range(len(rc))])
	#print "BL DEBUG:: r and imol is", r, imol

	self.failIf(r > 0.66, "Bad placement of terminal residue")

	# now test that the new atoms have the correct
	# B factor.
	new_atoms = residue_info(imol, "A", 0, "")
	self.failUnless(len(new_atoms) > 4,
                        "Not enough new atoms %s" %new_atoms)
	test_ls = map(lambda atom: atom[1][1], new_atoms)
	map(lambda atom: self.failUnlessAlmostEqual(atom[1][1], 45, 1,
                                                    "Fail b-factor test %s" %new_atoms),
            new_atoms)


    def test11_1(self):
        """Adding residue by phi psi, no crash"""

        imol = unittest_pdb("frag-2wot.pdb")
        self.failUnless(valid_model_molecule_qm(imol))
        v1 = add_terminal_residue_using_phi_psi(imol, "A", 275, "ALA", -60, -60)
        self.failUnlessEqual(v1, 1)
        v2 = add_terminal_residue_using_phi_psi(imol, "A", 276, "ALA", -60, -60)
        self.failUnlessEqual(v2, 1)
        v3 = add_terminal_residue_using_phi_psi(imol, "XX", 276, "ALA", -60, -60)
        self.failUnlessEqual(v3, 0)
        

    def test11_2(self):
        """Add Terminal Residue O Position"""

        imol = unittest_pdb("tutorial-modern.pdb")
        mtz_file_name = os.path.join(unittest_data_dir,
                                     "rnasa-1.8-all_refmac1.mtz")
        imol_map = make_and_draw_map(mtz_file_name, "FWT", "PHWT", "", 0, 0)

        # move the O close to where the N will end up - a bad place.
        # (we check that it moves there)
        #
        attribs = [[imol, "A", 93, "", " O  ", "", "x", 58.5],
                   [imol, "A", 93, "", " O  ", "", "y",  2.9],
                   [imol, "A", 93, "", " O  ", "", "z", -1.9]]
        set_atom_attributes(attribs)
        O_atom_o = get_atom(imol, "A", 93, "", " O  ", "")
        with NoBackups(imol):
            dummy = "dummy"
            add_terminal_residue(imol, "A", 93, "ALA", 1)
            N_atom_n = get_atom(imol, "A", 94, "", " N  ", "")
            O_atom_n = get_atom(imol, "A", 93, "", " O  ", "")
            dd_1 = bond_length_from_atoms(O_atom_o, N_atom_n)
            dd_2 = bond_length_from_atoms(O_atom_n, N_atom_n)
            print "Add terminal residue bond check dd_1", dd_1
            print "Add terminal residue bond check dd_2", dd_2
            
            # the new N will not always go into the same place
            # allow a bit more generosity in position difference
            self.failUnless(dd_1 < 0.4) # 0.4 should be generous enough
            self.failUnless(dd_2 > 1.9) # hooray

        
    def test12_0(self):
        """Select by Sphere"""

        global imol_rnase
        imol_sphere = new_molecule_by_sphere_selection(imol_rnase, 
                                                       24.6114959716797, 24.8355808258057, 7.43978214263916,
                                                       3.6, 1)

        self.failUnless(valid_model_molecule_qm(imol_sphere),
                        "Bad sphere molecule")
        n_atoms = 0
        for chain_id in chain_ids(imol_sphere):
            # add the number of atoms in the chains
            n_residues = chain_n_residues(chain_id, imol_sphere)
            print "   Sphere mol: there are %s residues in chain %s" %(n_residues, chain_id)
            for serial_number in range(n_residues):
                res_name = resname_from_serial_number(imol_sphere, chain_id, serial_number)
                res_no   = seqnum_from_serial_number(imol_sphere, chain_id, serial_number)
                ins_code = insertion_code_from_serial_number(imol_sphere, chain_id, serial_number)
                residue_atoms_info = residue_info(imol_sphere, chain_id, res_no, ins_code)
                n_atoms += len(residue_atoms_info)

        print "Found %s sphere atoms" %n_atoms
        self.failUnlessEqual(n_atoms, 20)


    def test13_0(self):
	"""Test Views"""

	view_number = add_view([32.0488, 21.6531, 13.7343],
			       [-0.12784, -0.491866, -0.702983, -0.497535],
			       20.3661,
			       "B11 View")
	go_to_view_number(view_number, 1)
	# test for something??


    def test13_1(self):
        """Delete Residue"""

        imol = read_pdb(rnase_pdb())
        self.failUnless(valid_model_molecule_qm(imol))
        delete_residue(imol, "A", 42, "")
        r = residue_info(imol, "A", 42, "")
        print "residue info (should be False):", r
        self.failIf(r)


    def test14_0(self):
	"""Label Atoms and Delete"""

	imol_frag = new_molecule_by_atom_selection(imol_rnase, "//B/10-12")
	set_rotation_centre(31.464, 21.413, 14.824)
	map(lambda n: label_all_atoms_in_residue(imol_frag, "B", n, ""), [10, 11, 12])
	rotate_y_scene(rotate_n_frames(200), 0.1)
	delete_residue(imol_frag, "B", 10, "")
	delete_residue(imol_frag, "B", 11, "")
	delete_residue(imol_frag, "B", 12, "")
	rotate_y_scene(rotate_n_frames(200), 0.1)
	# ???? what do we sheck for?


    def test15_0(self):
        """Rotamer outliers"""

        # pre-setup so that residues 1 and 2 are not rotamers "t" but 3
        # to 8 are (1 and 2 are more than 40 degrees away). 9-16 are
        # other residues and rotamers.
        imol_rotamers = read_pdb(os.path.join(unittest_data_dir,
                                              "rotamer-test-fragment.pdb"))

        rotamer_anal = rotamer_graphs(imol_rotamers)
        self.failUnless(type(rotamer_anal) is ListType)
        self.failUnlessEqual(len(rotamer_anal), 14)
        a_1 = rotamer_anal[0]
        a_2 = rotamer_anal[1]
        a_last = rotamer_anal[len(rotamer_anal)-1]

        anal_str_a1 = a_1[-1]
        anal_str_a2 = a_2[-1]
        anal_str_a3 = a_last[-1]

        # with new Molprobity rotamer probabilities,
        # residues 1 and 2 are no longer "not recognised"
        # they are in fact, just low probabilites.
        # 
        self.failUnless((anal_str_a1 == "VAL" and
                         anal_str_a2 == "VAL" and
                         anal_str_a3 == "Missing Atoms"),
                        "  failure rotamer test: %s %s %s" %(a_1, a_2, a_last))

        # Now we test that the probabilites of the
        # rotamer is correct:
        pr_1 = rotamer_anal[0][3]
        pr_2 = rotamer_anal[1][3]
        self.failUnless((pr_1 < 0.3 and pr_1 > 0.0),
                        "Failed rotamer outlier test for residue 1")
        self.failUnless((pr_2 < 0.3 and pr_2 > 0.0),
                        "Failed rotamer outlier test for residue 2")


    def test15_1(self):
        """HIS with unusual atom order rotates correct fragment for 180 sidechain flip"""

        imol = unittest_pdb("eleanor-HIS.pdb")
        self.failUnless(valid_model_molecule_qm(imol),
                        "BAD imol for 180 sidechain flip test imol: %s" %imol)

        N_atom_o = get_atom(imol, "A", 111, "", " N  ", "")
        ND1_Atom_o = get_atom(imol, "A", 111, "", " ND1", "")
        with NoBackups(imol):
            do_180_degree_side_chain_flip(imol, "A", 111, "", "")
            N_atom_n = get_atom(imol, "A", 111, "", " N  ", "")
            ND1_Atom_n = get_atom(imol, "A", 111, "", " ND1", "")

            # the N-atom stays still
            # the ND1 atom moves by > 1A.

            dd_1 = bond_length_from_atoms(N_atom_o, N_atom_n)
            dd_2 = bond_length_from_atoms(ND1_Atom_o, ND1_Atom_n)

            print "dd_1: %s dd_2: %s" %(dd_1, dd_2)

            self.failIf(dd_1 > 0.01,
                        "N atom moved - fail\n")

            self.failIf(dd_2 < 1.01,
                            "ND1 atom did not move enough - fail\n")


    # Don't reset the occupancies of the other parts of the residue
    # on regularization using alt confs
    #
    def test16_0(self):
        """Alt Conf Occ Sum Reset"""

        def get_occ_sum(imol_frag):

            def occ(att_atom_name, att_alt_conf, atom_ls):

                self.failUnless(type(atom_ls) is ListType) # repetition??
                for atom in atom_ls:
                    atom_name = atom[0][0]
                    alt_conf = atom[0][1]
                    occ = atom[1][0]
                    if ((atom_name == att_atom_name) and
                        (alt_conf == att_alt_conf)):
                        print "   For atom %s %s returning occ %s" %(att_atom_name, att_alt_conf, occ)
                        return occ

            # get_occ_sum body
            atom_ls = residue_info(imol_frag, "X", 15, "")
            self.failUnless(type(atom_ls) is ListType)
            return sum([occ(" CE ", "A", atom_ls), occ(" CE ", "B", atom_ls),
                        occ(" NZ ", "A", atom_ls), occ(" NZ ", "B", atom_ls)])

        # main body
        #
        imol_fragment = read_pdb(os.path.join(unittest_data_dir, "res098.pdb"))

        self.failUnless(valid_model_molecule_qm(imol_fragment),
                        "bad molecule for reading coords in Alt Conf Occ test")
        occ_sum_pre = get_occ_sum(imol_fragment)

        replace_state = refinement_immediate_replacement_state()
        set_refinement_immediate_replacement(1)
        regularize_zone(imol_fragment, "X", 15, 15, "A")
        accept_regularizement()
        if (replace_state == 0):
            set_refinement_immediate_replacement(0)

        occ_sum_post = get_occ_sum(imol_fragment)
        self.failUnlessAlmostEqual(occ_sum_pre, occ_sum_post, 1,
                                   "   test for closeness: %s %s" %(occ_sum_pre, occ_sum_post))


    def test16_1(self):
        """Correct occupancies after auto-fit rotamer on alt-confed residue"""

        imol = unittest_pdb("tutorial-modern.pdb")
        imol_map = make_and_draw_map(rnase_mtz(), "FWT", "PHWT", "", 0, 0)
        new_alt_conf = add_alt_conf(imol, "A", 93, "", "", 0)

        accept_regularizement()  # Presses the OK button for the alt conf

        auto_fit_best_rotamer(93, "A", "", "A", imol, imol_map, 1, 0.01)

        # Now test the ocupancies.  The problem was that we had been
        # ending up with atoms in the residue with occupancies of 0.0.
        # (They should be 0.8 and 0.2 - for me at least).  So test
        # that the atoms have occupancies of greater than 0.1 (I
        # suppose also test for occs > 0.85 would be good too).
        #
        atoms = residue_info(imol, "A", 93, "")
        for atom in atoms:
            occupancy = atom[1][0]
            alt_conf  = atom[0][1]
            self.failIf(occupancy < 0.1,
                        "bad occupancy in atom: %s" %atom)
            self.failIf(occupancy > 0.85,
                        "bad occupancy in atom: %s" %atom)


    def test16_2(self):
        """Rotamers work on MSE"""

        from types import ListType
        imol = unittest_pdb("pdb3knw.ent")
        se_1 = get_atom(imol, "A", 89, "", "SE  ")
        set_residue_to_rotamer_number(imol, "A", 89, "", "", 3)
        se_2 = get_atom(imol, "A", 89, "", "SE  ")

        print "    se_1:", se_1
        print "    se_2:", se_2

        self.failUnless(type(se_1) is ListType, "   could not find SE in pdb3knw.ent")
        self.failIf(atoms_match_qm(se_1, se_2))  # the SE moved, test passes
        
        
    def test16_3(self):
        """Hs are correctly swapped on a TYR"""

        imol = unittest_pdb("pdb1py3.ent")
        self.failUnless(valid_model_molecule_qm(imol),
                        "missing or bar pdb1py3")

        # the pdb file contains hydrogens and nomenclature
        # errors, lets see if we can fix them
        #
        fix_nomenclature_errors(imol)
        atoms = residue_info(imol, "C", 54, "")
        self.failUnless(isinstance(atoms, list), "atoms not a list")
        cd1 = get_atom_from_residue(" CD1", atoms, "")
        cd2 = get_atom_from_residue(" CD2", atoms, "")
        hd1 = get_atom_from_residue(" HD1", atoms, "")
        hd2 = get_atom_from_residue(" HD2", atoms, "")
        ce1 = get_atom_from_residue(" CE1", atoms, "")
        ce2 = get_atom_from_residue(" CE2", atoms, "")
        he1 = get_atom_from_residue(" HE1", atoms, "")
        he2 = get_atom_from_residue(" HE2", atoms, "")
        bonded_atoms = [[cd1, hd1],
                        [cd2, hd2],
                        [ce1, he1],
                        [ce2, he2]]
        results = map(lambda atom:
                      bond_length_within_tolerance_qm(atom[0], atom[1],
                                                      0.93, 0.02),
                      bonded_atoms)
        print "BL DEBUG:: results:", results
        self.failUnless(all(results))


    def test16_4(self):
        """Splitting residue leaves no atoms with negative occupancy"""

        # return False if there are negative occupancies
        #
        def check_for_negative_occs(occs):
            if not occs:    # check for list?
                return False
            for occ in occs:
                if occ < 0.0:
                    return False
            return True

        imol = unittest_pdb("tutorial-modern.pdb")
        mtz_file_name = rnase_mtz()
        imol_map = make_and_draw_map(mtz_file_name, "FWT", "PHWT", "", 0, 0)

        zero_occupancy_residue_range(imol, "A", 37, 37)
        new_alt_conf = add_alt_conf(imol, "A", 37, "", "", 0)
        atoms = residue_info(imol, "A", 37, "")
        occs = [atom[1][0] for atom in atoms]
        if (len(occs) < 5):
            return False # too few atoms
        else:
            occs_ok_status = check_for_negative_occs(occs)
            self.failUnless(occs_ok_status,
                            "Ooops: bad occupancies: %s" %occs)
            close_molecule(imol)
            close_molecule(imol_map)


    def test17_0(self):
        """Pepflip flips the correct alt confed atoms"""

        imol = unittest_pdb("alt-conf-pepflip-test.pdb")
        self.failUnless(valid_model_molecule_qm(imol))

        # get the originla coords
        c_atom_A_o = get_atom(imol, "A", 65, "", " C  ", "A")
        o_atom_A_o = get_atom(imol, "A", 65, "", " O  ", "A")
        n_atom_A_o = get_atom(imol, "A", 66, "", " N  ", "A")
        c_atom_B_o = get_atom(imol, "A", 65, "", " C  ", "B")
        o_atom_B_o = get_atom(imol, "A", 65, "", " O  ", "B")
        n_atom_B_o = get_atom(imol, "A", 66, "", " N  ", "B")

        pepflip(imol, "A", 65, "", "B")

        # get the new coords
        c_atom_A_n = get_atom(imol, "A", 65, "", " C  ", "A")
        o_atom_A_n = get_atom(imol, "A", 65, "", " O  ", "A")
        n_atom_A_n = get_atom(imol, "A", 66, "", " N  ", "A")
        c_atom_B_n = get_atom(imol, "A", 65, "", " C  ", "B")
        o_atom_B_n = get_atom(imol, "A", 65, "", " O  ", "B")
        n_atom_B_n = get_atom(imol, "A", 66, "", " N  ", "B")

        # now, the *A-n atoms should match the position of the
        # *A-o atoms:
        b1 = bond_length_from_atoms(c_atom_A_o, c_atom_A_n)
        b2 = bond_length_from_atoms(o_atom_A_o, o_atom_A_n)
        b3 = bond_length_from_atoms(n_atom_A_o, n_atom_A_n)
        b4 = bond_length_from_atoms(c_atom_B_o, c_atom_B_n)
        b5 = bond_length_from_atoms(o_atom_B_o, o_atom_B_n)
        b6 = bond_length_from_atoms(n_atom_B_o, n_atom_B_n)

        self.failUnless(b1 < 0.001 and b2 < 0.001 and b3 < 0.001,
                        "   bad! A conf moved %s %s %s" %(b1, b2, b3))

        self.failUnless(b4 > 0.8 and b5 > 2.0 and b6 > 0.5,
                        "   bad! B conf moved too little %s %s %s" %(b4, b5, b6))


    # This test we expect to fail until the CISPEP correction code is in
    # place (using mmdb-1.10+).
    # 
    def test18_0(self):
        """Correction of CISPEP test"""

        # In this test the cis-pep-12A has indeed a CIS pep and has been
        # refined with refmac and is correctly annotated.  There was 4 on
        # read.  There should be 3 after we fix it and write it out.
        # 
        # Here we sythesise user actions:
        # 1) CIS->TRANS on the residue 
        # 2) convert leading residue to a ALA
        # 3) mutate the leading residue and auto-fit it to correct the chiral
        #    volume error :)

        chain_id = "A"
        resno = 11
        ins_code = ""

        cis_pep_mol = read_pdb(os.path.join(unittest_data_dir,
                                            "tutorial-modern-cis-pep-12A_refmac0.pdb"))
        if self.skip_test(cis_pep_mol < 0, 
                          "skipping CIS test as problem on read pdb"):
            return

        view_number = add_view([63.455, 11.764, 1.268],
                               [-0.760536, -0.0910907, 0.118259, 0.631906],
                               15.7374,
                               "CIS-TRANS cispep View")
        go_to_view_number(view_number, 1)

        cis_trans_convert(cis_pep_mol, chain_id, resno, ins_code)
        pepflip(cis_pep_mol, chain_id, resno, ins_code, "")
        res_type = residue_name(cis_pep_mol, chain_id, resno, ins_code)
        self.failUnless(type(res_type) is StringType)
        mutate(cis_pep_mol, chain_id, resno, "", "GLY")
        with_auto_accept(
                [refine_zone, cis_pep_mol, chain_id, resno, (resno + 1), ""],
                #[accept_regularizement],
                [mutate, cis_pep_mol, chain_id, resno, "", res_type],
                [auto_fit_best_rotamer, resno, "", ins_code, chain_id, cis_pep_mol,
                 imol_refinement_map(), 1, 1],
        [with_auto_accept, [refine_zone, cis_pep_mol, chain_id, resno, (resno + 1), ""]]
                #[accept_regularizement]
                )

        # replace at some point with call_with_input_file (once we have that)
        # actually already there using "with file open"...
        tmp_file = "tmp-fixed-cis.pdb"
        write_pdb_file(cis_pep_mol, tmp_file)
        # assuming we do not always have grep we extract information with
        # a python funcn
        # first arg: search string
        # second arg: filename
        def grep_to_list(pattern_str, filen):
            import re
            ret = []
            go = False
            try:
                fin = open(filen, 'r')
                lines = fin.readlines()
                fin.close()
                go = True
            except IOError:
                print "file no found"
            if go:
                pattern = re.compile(pattern_str)
                for line in lines:
                    match = pattern.search(line)
                    if match:
                        ret.append(line)
            return ret

        o = grep_to_list("CISPEP", tmp_file)
        self.failUnlessEqual(len(o), 3)


    def test18_1(self):
        """H on a N moves on cis-trans convert"""

        imol = unittest_pdb("tutorial-modern.pdb")
        imol_2 = new_molecule_by_atom_selection(imol, "//A/1-10")
        coot_reduce(imol_2)
        H_atom_o = get_atom(imol_2, "A", 6, "", " H  ", "")
        with NoBackups(imol_2):
            cis_trans_convert(imol_2, "A", 5, "") # 5-6 peptide
        H_atom_n = get_atom(imol_2, "A", 6, "", " H  ", "")
        dd = bond_length_from_atoms(H_atom_o, H_atom_n)
        close_molecule(imol)
        close_molecule(imol_2)
        print "dd:", dd
        self.failUnless(dd> 1.4)


    def test18_2(self):
        """HA on a ALA exists after mutation to GLY"""

        imol = unittest_pdb("tutorial-modern.pdb")
        self.failUnless(valid_model_molecule_qm(imol))
        imol_2 = new_molecule_by_atom_selection(imol, "//A/5-11")
        self.failUnless(valid_model_molecule_qm(imol_2))
        coot_reduce(imol_2)
        H_atom_o = get_atom(imol_2, "A", 10, "", " HA ", "")
        self.failUnless(isinstance(H_atom_o, list))
        mutate(imol_2, "A", 10, "", "GLY")
        H_atom_n = get_atom(imol_2, "A", 10, "", " HA ", "")
        self.failIf(isinstance(H_atom_n, list),
                    "atom still exists %s" %H_atom_n)


    def test19_0(self):
        """Refine Zone with Alt conf"""

        imol = unittest_pdb("tutorial-modern.pdb")
        mtz_file_name = (os.path.join(unittest_data_dir, "rnasa-1.8-all_refmac1.mtz"))
        imol_map = make_and_draw_map(mtz_file_name, "FWT", "PHWT", "", 0, 0)

        set_imol_refinement_map(imol_map)
        at_1 = get_atom(imol, "B", 72, "", " SG ", "B")
        with_auto_accept([refine_zone, imol, "B", 72, 72, "B"])
        at_2 = get_atom(imol, "B", 72, "", " SG ", "B")
        d = bond_length_from_atoms(at_1, at_2)
        # the atom should move in the refinement
        print "   refined moved: d=", d
        # 20120110 new-style NBCs means that the
        # atoms move less here
	# 20160608 - they move still less  (not sure why this time)
	#
        self.failIf(d < 0.09,
                    "   refined atom failed to move: d=%s" %d)


    def test20_0(self):
        """Sphere Refine"""

        def sphere_refine_here():
            from types import ListType
            active_atom = active_residue()
            self.failUnless(type(active_atom) is ListType)
            centred_residue = active_atom[1:4]
            imol = active_atom[0]
            other_residues = residues_near_residue(imol,
                                                   centred_residue,
                                                   3.2)
            all_residues = [centred_residue]
            if other_residues:
                all_residues += other_residues
            print "   imol: %s residues: %s" %(imol, all_residues)
            with_auto_accept([refine_residues, imol, all_residues])

        # main line
        imol = unittest_pdb("tutorial-add-terminal-1-test.pdb")
        self.failUnless(valid_model_molecule_qm(imol),
                        "   molecule pdb not found")

        new_alt_conf = add_alt_conf(imol, "A", 93, "", "", 0)

        set_go_to_atom_molecule(imol)
        success = set_go_to_atom_chain_residue_atom_name("A", 93, " CA ")
        self.failIf(success == 0, "   failed to go to A93 CA atom")
        sphere_refine_here()

        set_go_to_atom_chain_residue_atom_name("A", 41, " CA ")
        sphere_refine_here()

        # these should be a standard peptide bond length - however, if we don't find
        # the peptide link, they get pushed apart.  Test for that.
        atom_1 = get_atom(imol, "A", 40, "", " C  ", "")
        atom_2 = get_atom(imol, "A", 41, "", " N  ", "")
        bl = bond_length_from_atoms(atom_1, atom_2)
        print "======= got bond length", bl
        self.failUnless(bl < 1.4)


    def test20_1(self):
        """Refinement gives useful results"""

#        if self.skip_test(True, "skipping because we dont get useful refinement " \
#                          "results due to threading!"):
#            return
        #
        def no_non_bonded(ls):
            ret = []
            for item in ls:
                if not (item[0] == "Non-bonded"):
                    ret.append(item)
            return ret

        # return False or a number, which is the current overweighting of
        # the density terms.  When not overweighted, the geometric chi
        # squareds will be 1.0.
        #
        def weight_scale_from_refinement_results(rr):
            self.failUnless(rr, "refinement returned False")
            nnb_list = no_non_bonded(rr[2])
            chi_squares = [x[2] for x in nnb_list]
            n = len(chi_squares)
            summ = sum(chi_squares)
            return summ / n

        imol = read_pdb(rnase_pdb())
        imol_map = make_and_draw_map(rnase_mtz(), "FWT", "PHWT", "", 0, 0)
        set_imol_refinement_map(imol_map)

        # not convinced this should be an 'indefinite' loop,
        # so we just try 10 times (wild guess for now)
        failed = True
        for idum in range(10):
            refine_zone_with_full_residue_spec(imol, "A", 40, "", 43, "", "")
            results = accept_moving_atoms()
            print "   refinement results:", results
            ow = weight_scale_from_refinement_results(results)
            print "::::   ow factor", ow
            if (ow < 1.1 and ow > 0.9):
                failed = False
                break
            else:
                new_weight = matrix_state() / (ow * ow)
                print "   INFO:: setting refinement weight to", new_weight, "from", matrix_state() , "/ (", ow, '*', ow, ')'
                set_matrix(new_weight)

        # this test doesn't converge on Ubuntu for some reason that I don't understand
        failed = False

        self.failIf(failed,
                    "Refinement didnt 'converge' in 5 rounds")


    def test20_2(self):
        """Neighbour-Refine doesn't destroy disulfide bonds"""

        global imol_rnase_map
        
        imol = unittest_pdb("tutorial-modern.pdb")

        set_imol_refinement_map(imol_rnase_map)
        set_go_to_atom_molecule(imol)
        set_go_to_atom_chain_residue_atom_name("B", 7, " CA ")

        rc_spec = ["B", 7, ""]
        ls = residues_near_residue(imol, rc_spec, 2.2)
        residue_list = rc_spec + ls

        with AutoAccept():
            refine_residues(imol, residue_list)
        at_spec_1 = ["B",  7, "", " SG ", ""]
        at_spec_2 = ["B",  96, "", " SG ", ""]
        at_1 = get_atom_from_spec(imol, at_spec_1)
        at_2 = get_atom_from_spec(imol, at_spec_2)
        bl = bond_length_from_atoms(at_1, at_2)

        state = bond_length_within_tolerance_qm(at_1, at_2, 2.0, 0.05)

        self.failUnless(state)

        # do it again
        with AutoAccept():
            refine_residues(imol, residue_list)
        at_spec_1 = ["B",  7, "", " SG ", ""]
        at_spec_2 = ["B",  96, "", " SG ", ""]
        at_1 = get_atom_from_spec(imol, at_spec_1)
        at_2 = get_atom_from_spec(imol, at_spec_2)
        bl = bond_length_from_atoms(at_1, at_2)

        state = bond_length_within_tolerance_qm(at_1, at_2, 2.0, 0.05)

        self.failUnless(state)        

        
    def test21_0(self):
        """Rigid Body Refine Alt Conf Waters"""

        imol_alt_conf_waters = read_pdb(os.path.join(unittest_data_dir, "alt-conf-waters.pdb"))

        rep_state = refinement_immediate_replacement_state()
        set_refinement_immediate_replacement(1)
        refine_zone(imol_alt_conf_waters, "D", 71, 71, "A")
        accept_regularizement()
        refine_zone(imol_alt_conf_waters, "D", 2, 2, "")
        accept_regularizement()
        set_refinement_immediate_replacement(rep_state)
        # what do we check for?!


    def test22_0(self):
        """Setting multiple atom attributes"""

        global imol_rnase
        self.failUnless(valid_model_molecule_qm(imol_rnase),
                        "   Error invalid imol-rnase")

        x_test_val = 2.1
        y_test_val = 2.2
        z_test_val = 2.3
        o_test_val = 0.5
        b_test_val = 44.4
        ls = [[imol_rnase, "A", 2, "", " CA ", "", "x", x_test_val],
              [imol_rnase, "A", 2, "", " CA ", "", "y", y_test_val],
              [imol_rnase, "A", 2, "", " CA ", "", "z", z_test_val],
              [imol_rnase, "A", 2, "", " CA ", "", "occ", o_test_val],
              [imol_rnase, "A", 2, "", " CA ", "", "b", b_test_val]]

        set_atom_attributes(ls)
        atom_ls = residue_info(imol_rnase, "A", 2, "")
        self.failIf(atom_ls == [])
        for atom in atom_ls:
            atom_name = atom[0][0]
            if (atom_name == " CA "):
                x = atom[2][0]
                y = atom[2][1]
                z = atom[2][2]
                occ = atom[1][0]
                b = atom[1][1]
                self.failUnlessAlmostEqual(x, x_test_val, 1,
                                           "Error in setting multiple atom: These are not close %s %s" %(x, x_test_val))
                self.failUnlessAlmostEqual(y, y_test_val, 1,
                                           "Error in setting multiple atom: These are not close %s %s" %(y, y_test_val))
                self.failUnlessAlmostEqual(z, z_test_val, 1,
                                           "Error in setting multiple atom: These are not close %s %s" %(z, z_test_val))
                self.failUnlessAlmostEqual(occ, o_test_val, 1,
                                           "Error in setting multiple atom: These are not close %s %s" %(occ, o_test_val))
                self.failUnlessAlmostEqual(b, b_test_val, 1,
                                           "Error in setting multiple atom: These are not close %s %s" %(b, b_test_val))


    def test23_0(self):
        """Tweak Alt Confs on Active Residue"""

        global imol_rnase
        # did it reset?
        def matches_alt_conf_qm(imol, chain_id, resno, inscode, atom_name_ref, alt_conf_ref):
            atom_ls = residue_info(imol, chain_id, resno, inscode)
            self.failIf(atom_ls ==[],
                        "No atom list found - failing.")
            for atom in atom_ls:
                atom_name = atom[0][0]
                alt_conf = atom[0][1]
                if ((atom_name == atom_name_ref) and (alt_conf == alt_conf_ref)):
                    return True
            return False

        # main line
        chain_id = "B"
        resno = 58
        inscode = ""
        atom_name = " CB "
        new_alt_conf_id = "B"
        set_go_to_atom_molecule(imol_rnase)
        set_go_to_atom_chain_residue_atom_name(chain_id, resno, " CA ")
        set_atom_string_attribute(imol_rnase, chain_id, resno, inscode,
                                  atom_name, "", "alt-conf", new_alt_conf_id)
        self.failUnless(matches_alt_conf_qm(imol_rnase, chain_id, resno, inscode,
                                            atom_name, new_alt_conf_id),
                        "   No matching pre CB altconfed - failing.")
        sanitise_alt_confs_in_residue(imol_rnase, chain_id, resno, inscode)
        self.failUnless(matches_alt_conf_qm(imol_rnase, chain_id, resno, inscode, atom_name, ""),
                        "   No matching post CB (unaltconfed) - failing.")


    def test24_0(self):
        """Backrub rotamer"""

        global imol_rnase_map
        #
        def get_mover_list(imol_1, imol_2, res_no):
            atoms_1 = residue_info(imol_1, "A", res_no, "")
            atoms_2 = residue_info(imol_2, "A", res_no, "")
            mover_list = []
            for atom_1, atom_2 in zip(atoms_1, atoms_2):
                pos_1 = atom_1[2]
                pos_2 = atom_2[2]
                ls = map(close_float_qm, pos_1, pos_2)
                if not all(ls):
                    # did move
                    mover_list.append(atom_1[0][0])
            mover_list.sort()
            return mover_list

        # main line
        set_imol_refinement_map(imol_rnase_map)

        imol = handle_read_draw_molecule_with_recentre(os.path.join(unittest_data_dir,
                                                                    "backrub-fragment.pdb"),
                                                       0)
        imol_copy = copy_molecule(imol)

        self.failUnless(valid_model_molecule_qm(imol),
                        "backrub-rotamer backrub-fragment.pdb not found it seems")

        status = backrub_rotamer(imol_copy, "A", 37, "", "")

        self.failUnless(status == 1,
                        "backrub-rotamer status failure")

        # Now, N should have moved in 38 (nothing else)
        #
        # C and O should have moved in 36 (nothing else)
        # 
        # all atoms in 37.
        # 

        atoms_37 = residue_info(imol, "A", 37, "")
        calc_37_movers = [atom[0][0] for atom in atoms_37]
        calc_37_movers.sort()
        calc_38_movers = [" N  "]
        calc_36_movers = [" C  ", " O  "]

        movers_36 = get_mover_list(imol, imol_copy, 36)
        movers_37 = get_mover_list(imol, imol_copy, 37)
        movers_38 = get_mover_list(imol, imol_copy, 38)

        #print "calc_36_movers:", calc_36_movers
        #print " obs_36_movers:",         movers_36
        #print "calc_37_movers:", calc_37_movers
        #print " obs_37_movers:",         movers_37
        #print "calc_38_movers:", calc_38_movers
        #print " obs_38_movers:",         movers_38

        self.failUnlessEqual(movers_36, calc_36_movers,
                             "  fail on 36 movers similar")

        self.failUnlessEqual(movers_37, calc_37_movers,
                             "  fail on 37 movers similar")

        self.failUnlessEqual(movers_38, calc_38_movers,
                             "  fail on 36 movers similar")


    def test25_0(self):
        """Libcif horne"""

        if self.skip_test(not os.path.isfile(horne_pdb),
                          "file %s not found - skipping test" %horne_pdb):
            return
        if self.skip_test(not os.path.isfile(horne_cif),
                          "file %s not found - skipping test" %horne_cif):
            return
        if self.skip_test(not os.path.isfile(horne_works_cif),
                          "file %s not found - skipping test" %horne_works_cif):
            return

        imol = read_pdb(horne_pdb)
        self.failUnless(valid_model_molecule_qm(imol),
                        "bad molecule number %i" %imol)

        read_cif_dictionary(horne_works_cif)
        with_auto_accept(
                [regularize_zone, imol, "A", 1, 1, ""]
                )
        print "\n \n \n \n"
        read_cif_dictionary(horne_cif)
        with_auto_accept(
                [regularize_zone, imol, "A", 1, 1, ""]
                )
        cent = molecule_centre(imol)
        self.failUnless(cent[0] < 2000,
                        "Position fail 2000 test: %s in %s" %(cent[0], cent))


    def test26_0(self):
        """Refmac Parameters Storage"""

        arg_list = [rnase_mtz(), "/RNASE3GMP/COMPLEX/FWT", "/RNASE3GMP/COMPLEX/PHWT", "", 0, 0, 1, 
                    "/RNASE3GMP/COMPLEX/FGMP18", "/RNASE3GMP/COMPLEX/SIGFGMP18",
                    "/RNASE/NATIVE/FreeR_flag", 1]
        imol = make_and_draw_map_with_refmac_params(*arg_list)

        self.failUnless(valid_map_molecule_qm(imol),
                        "  Can't get valid refmac map molecule from %s" %rnase_mtz())

        refmac_params = refmac_parameters(imol)
        refmac_params[0] = os.path.normpath(refmac_params[0])

        self.failUnlessEqual(refmac_params, arg_list,
                             "        non matching refmac params")


    def test26_1(self):
        """OXT is added before TER record - add only one"""

        imol = unittest_pdb("test-TER-OXT.pdb")
        opdb = "test-TER-OXT-added.pdb"

        self.failUnless(valid_model_molecule_qm(imol),
                        "Failed to read test-TER-OXT.pdb")

        add_status_1 = add_OXT_to_residue(imol, "A", 14, "")
        add_status_2 = add_OXT_to_residue(imol, "A", 14, "")

        # the second add should be ignored, only 1 OXT allowed
        # per residue.

        self.failUnless(add_status_1 == 1,
                        "   add-status-1 not success - fail")

        self.failUnless(add_status_2 == 0,
                        "   add-status-2 was success - fail")

        write_pdb_file(imol, opdb)
        self.failUnless(os.path.isfile(opdb),
                        "Failed to find test-TER-OXT-added.pdb")
        fin = open(opdb, 'r')
        lines = fin.readlines()
        fin.close()
        count = 0
        ter_line = False
        oxt_line = False
        for line in lines:
            self.failIf(line == lines[-1],
                        "EOF and no TER and OXT")
            self.failIf(line[0:3] == "END",
                        "Found END before TER and OXT")
            if (line[0:3] == "TER"):
                print "   found TER", line
                ter_line = True
                self.failUnless(oxt_line)   # TER happens after OXT line
                break # pass (have TER after OXT)
            if (line[13:16] == "OXT"):
                self.failIf(ter_line)  # fail because TER has already happened
                self.failIf(oxt_line, "   Encountered another OXT! - fail")
                oxt_line = True


    def test27_0(self):
        """The position of the oxygen after a mutation"""

        import operator
        # Return the o-pos (can be False) and the number of atoms read.
        #
        def get_o_pos(pdb_file):
            fin = open(pdb_file, 'r')
            lines = fin.readlines()
            fin.close
            n_atoms = 0
            o_pos = False
            for line in lines:
                atom_name = line[12:16]
                if (line[0:4] == "ATOM"):
                    n_atoms +=1
                if (atom_name == " O  "):
                    o_pos = n_atoms

            return o_pos, n_atoms

        # main body
        #
        hist_pdb_name = "his.pdb"

        imol = read_pdb(os.path.join(unittest_data_dir, "val.pdb"))
        self.failUnless(valid_model_molecule_qm(imol),
                        "   failed to read file val.pdb")
        mutate_state = mutate(imol, "C", 3, "", "HIS")
        self.failUnlessEqual(mutate_state, 1,
                             "   failure in mutate function")
        write_pdb_file(imol, hist_pdb_name)
        self.failUnless(os.path.isfile(hist_pdb_name),
                        "   file not found: %s" %hist_pdb_name)

        o_pos, n_atoms = get_o_pos(hist_pdb_name)
        self.failUnlessEqual(n_atoms, 10,
                             "   found %s atoms (not 10)" %n_atoms)
        self.failUnless(operator.isNumberType(o_pos),
                        "   Oxygen O position not found")
        self.failUnlessEqual(o_pos, 4,
                             "   found O atom at %s (not 4)" %o_pos) 


    def test27_1(self):
        """TER is at the end of a nucleotide after mutation"""

        # Before this test, if we mutated a nucleotide at the end of a
        # chain, then the TER record appeared in the PDB file before the
        # additional new atoms.  Wrongness.  James Parker bug.

        imol = unittest_pdb("2yie-frag.pdb")
        self.failUnless(valid_model_molecule_qm(imol),
                        "   failed to read 2yie-frag.pdb")

        status = mutate_base(imol, "X", 54, "", "C")

        self.failUnlessEqual(status, 1,
                             "failed to mutate 2yie-frag.pdb")

        write_pdb_file(imol, "2yie-mutated.pdb")

        fin = open("2yie-mutated.pdb", 'r')
        lines = fin.readlines()
        fin.close()
        ter_line = False
        # check for END in last line
        self.failUnlessEqual(lines[-1][0:3], "END",
                             "'END' is not in the end")
        for line in lines:
            if (line[0:3] == "TER"):
                ter_line = True
            if (line[0:4] == "ATOM"):
                self.failIf(ter_line, "TER before ATOM")  # fail because TER has already happened


    def test27_2(self):
        """C7 is removed on mutation from a DC"""

        imol = unittest_pdb("4f8g.pdb")
        self.failUnless(valid_model_molecule_qm(imol),
                        "Missing 4fg8.pdb")

        status = mutate_base(imol, "A", 19, "", "DC")
        self.failUnlessEqual(status, 1,
                             "failed to mutate 1")

        atom_list = residue_info(imol, "A", 19, "")
        self.failUnless(isinstance(atom_list, list),
                        "Bad atom list 1")

        self.failUnless(len(atom_list) > 10,
                        "Bad atom list, need more atoms 1")

        atoms = map(residue_atom2atom_name, atom_list)
        self.failIf(" C7 " in atoms,
                    "C7 still present")

        write_pdb_file(imol, "test-1.pdb")  # why?


    def test27_3(self):
        """C7 is added on mutation to a DC"""

        imol = unittest_pdb("4f8g.pdb")
        # now try mutate 10 to a T

        status = mutate_base(imol, "A", 10, "" , "DT")
        self.failUnlessEqual(status, 1,
                             "failed to mutate back 4f8g.pdb")
        atom_list = residue_info(imol, "A", 10, "")
        self.failUnless(isinstance(atom_list, list),
                        "Bad atom list 2")

        atoms = map(residue_atom2atom_name, atom_list)
        self.failUnless(" C7 " in atoms,
                        "C7 not present in T")

        write_pdb_file(imol, "test-2.pdb")  # why?
        
    
    def test28_0(self):
        """Deleting (non-existing) Alt conf and Go To Atom [JED]"""

        global imol_rnase
        # alt conf "A" does not exist in this residue:
        delete_residue_with_full_spec(imol_rnase, 1, "A", 88, "", "A")
        # to activate the bug, we need to search over all atoms
        active_residue()
        # test for what?? (no crash??)


    def test29_0(self):
        """Mask and difference map"""

        def get_ca_coords(imol_model, resno_1, resno_2):
            coords = []
            for resno in range(resno_1, resno_2+1):
                atom_info = atom_specs(imol_model, "A", resno, "", " CA ", "")
                co = atom_info[3:6]
                coords.append(co)
            return coords

        d_1 = difference_map(-1, 2, -9)
        d_2 = difference_map(2, -1, -9)

        self.failUnlessEqual(d_1, -1, "failed on bad d_1")
        self.failUnlessEqual(d_2, -1, "failed on bad d_1")

        imol_map_nrml_res = make_and_draw_map(rnase_mtz(), "FWT", "PHWT", "", 0, 0)
        prev_sampling_rate = get_map_sampling_rate()
        nov_1 = set_map_sampling_rate(2.2)
        imol_map_high_res = make_and_draw_map(rnase_mtz(), "FWT", "PHWT", "", 0, 0)
        nov_2 = set_map_sampling_rate(prev_sampling_rate)
        imol_model = handle_read_draw_molecule_with_recentre(rnase_pdb(), 0)

        imol_masked = mask_map_by_atom_selection(imol_map_nrml_res, imol_model,
                                                 "//A/1-10", 0)

        self.failUnless(valid_map_molecule_qm(imol_map_nrml_res))

        # check that imol-map-nrml-res is good here by selecting
        # regions where the density should be positive. The masked
        # region should be 0.0
        #
        high_pts = get_ca_coords(imol_model, 20, 30)
        low_pts  = get_ca_coords(imol_model,  1, 10)

        high_values = map(lambda pt: density_at_point(imol_masked, *pt),
                          high_pts)
        low_values = map(lambda pt: density_at_point(imol_masked, *pt),
                         low_pts)

        print "   high-values: %s  low values: %s" %(high_values, low_values)

        self.failUnless((sum(low_values) < 0.000001), "Bad low values")

        self.failUnless((sum(high_values) > 5), "Bad high values")

        diff_map = difference_map(imol_masked, imol_map_high_res, 1.0)
        self.failUnless(valid_map_molecule_qm(diff_map),
                        "failure to make good difference map")

        # now in diff-map low pt should have high values and high
        # pts should have low values

        diff_high_values = map(lambda pt: density_at_point(diff_map, *pt),
                               high_pts)
        diff_low_values = map(lambda pt: density_at_point(diff_map, *pt),
                              low_pts)

        print "   diff-high-values: %s  diff-low-values: %s" %(diff_high_values, diff_low_values)

        sum_diff_high_values = sum(diff_high_values) 
        self.failUnless((sum_diff_high_values < 0.06),
                        "Bad diff high values: value %f target: %f"
                        %(sum_diff_high_values, 0.06))

        sum_diff_low_values = sum(diff_low_values)
        self.failUnless((sum_diff_low_values < -5),
                        "Bad diff low values %s" %sum_diff_low_values)


    def test29_1(self):
        """Skeletonize a map"""

        imol = read_pdb(rnase_pdb())
        imol_map = make_and_draw_map(rnase_mtz(), "FWT", "PHWT", "", 0, 0)
        # should check that both are ok.
        self.failUnless(valid_model_molecule_qm(imol),"failed to get imol")
        self.failUnless(valid_map_molecule_qm(imol_map), "failed to get imap")

        skeletonize_map(1, imol_map)
        skeletonize_map(0, imol_map)
        skeletonize_map(-1, -1)
        skeletonize_map(0, 0)
        close_molecule(imol)
        close_molecule(imol_map)
        # yeah, no crash
        

    def test29_2(self):
        """Simple Averaged maps"""

        imol_map_1 = make_and_draw_map(rnase_mtz(), "FWT", "PHWT", "", 0, 0)
        novalue_1  = set_map_sampling_rate(2.5)
        imol_map_2 = make_and_draw_map(rnase_mtz(), "FWT", "PHWT", "", 0, 0)
        novalue_2  = set_map_sampling_rate(1.5) # reset it

        new_map = average_map([[imol_map_1, 1],
                               [imol_map_2, 1]])

        self.failUnless(valid_map_molecule_qm(new_map),
                        " average map fail")

        set_ignore_pseudo_zeros_for_map_stats(0)

        # the differnce map shoudl be nearly flat (0.0)
        #
        diff_map = difference_map(imol_map_1, new_map, 1.0)

        rms_1 = map_sigma(new_map)
        rms_2 = map_sigma(diff_map)

        print "  INFO:: map sigmas: normal %s and diff-map: %s" %(rms_1, rms_2)
        self.failUnless(rms_1 > 0.3,
                        " map sigma 1 fail")

        self.failUnless(rms_2 < 0.003,
                        " map sigma for diff average map fail")


    def test30_0(self):
        """Make a glycosidic linkage"""

        carbo = "multi-carbo-coot-3.pdb"
        imol = unittest_pdb(carbo)

        if self.skip_test(not valid_model_molecule_qm(imol),
                          "file %s not found, skipping test" %carbo):
            return

        atom_1 = get_atom(imol, "A", 1, "", " O4 ")
        atom_2 = get_atom(imol, "A", 2, "", " C1 ")

        print "   bond-length: ", bond_length(atom_1[2], atom_2[2])

        s = dragged_refinement_steps_per_frame()
        set_dragged_refinement_steps_per_frame(300)
        with_auto_accept([regularize_zone, imol, "A", 1, 2, ""])
        set_dragged_refinement_steps_per_frame(s)

        atom_1 = get_atom(imol, "A", 1, "", " O4 ")
        atom_2 = get_atom(imol, "A", 2, "", " C1 ")

        print "   bond-length: ", bond_length(atom_1[2], atom_2[2])

    def test30_1(self):
        """Refine an NAG-ASN Link"""

        imol = unittest_pdb("pdb2qc1-sans-cho.pdb")

        if self.skip_test(not valid_model_molecule_qm(imol),
            "failed to find pdb2qc1-sans-cho.pdb, skipping test" ):
            return

        imol_map = make_and_draw_map(rnase_mtz(), "FWT", "PHWT", "", 0, 0)
        status = add_linked_residue(imol, "B", 141, "", "NAG", "NAG-ASN", 3000)
        # do something with status?!
        with_auto_accept([refine_residues, imol, [["B", 141, ""], ["B", 464, ""]]])

        atom_1 = get_atom(imol, "B", 141, "", " ND2")
        atom_2 = get_atom(imol, "B", 464, "", " C1 ")

        self.failUnless(bond_length_within_tolerance_qm(atom_1, atom_2, 1.43, 0.2),
                        "the bond between the two atoms is not within the tolerance")


    def test31_0(self):
        """Test for flying hydrogens on undo"""

        from types import ListType

        imol = unittest_pdb("monomer-VAL.pdb")

        self.failUnless(valid_model_molecule_qm(imol),
                        "   Failure to read monomer-VAL.pdb")

        with_auto_accept([regularize_zone, imol, "A", 1, 1, ""])
        set_undo_molecule(imol)
        apply_undo()
        with_auto_accept([regularize_zone, imol, "A", 1, 1, ""])

        atom_1 = get_atom(imol, "A", 1, "", "HG11")
        atom_2 = get_atom(imol, "A", 1, "", " CG1")

        self.failUnless((type(atom_1) is ListType) and
                        (type(atom_2) is ListType) ,
                        "   flying hydrogen failure, atoms: %s %s" %(atom_1, atom_2))
        self.failUnless(bond_length_within_tolerance_qm(atom_1, atom_2, 0.96, 0.02),
                        "   flying hydrogen failure, bond length %s, should be 0.96" %bond_length_from_atoms(atom_1, atom_2))


    def test32_0(self):
        """Test for regularization and mangling of hydrogen names from a PDB v 3.0"""

        # Note that it seems to me that the bonds are not within
        # tolerance before the regularization.
        # 
        # Also note that the D atom has been removed because now we have
        # a test for the atom names matching (hydrogens are not checked
        # currently, but this atom does not appear to be a hydrogen
        # (element is " D").

        imol = unittest_pdb("3ins-6B-3.0-no-peptide-D.pdb")
        self.failUnless(valid_model_molecule_qm(imol),
                        "Bad read of unittest test pdb: 3ins-6B-3.0-no-peptide-D.pdb")

        with_auto_accept([regularize_zone, imol, "B", 6, 6, ""])
        atom_pairs = [["HD11", " CD1"],
                      ["HD12", " CD1"],
                      ["HD13", " CD1"],
                      ["HD21", " CD2"],
                      ["HD22", " CD2"],
                      ["HD23", " CD2"]]
        atoms_1 = map(lambda pair: get_atom(imol, "B", 6, "", pair[0]), atom_pairs)
        atoms_2 = map(lambda pair: get_atom(imol, "B", 6, "", pair[1]), atom_pairs)
        #all_true = map(lambda atom_1, atom_2: bond_length_within_tolerance_qm(atom_1, atom_2, 0.96, 0.02), atoms_1, atoms_2)
        all_false = filter(lambda (atom_1, atom_2): not bond_length_within_tolerance_qm(atom_1, atom_2, 0.96, 0.02), zip(atoms_1, atoms_2))
        #all_true = filter(lambda (atom_1, atom_2): bond_length_within_tolerance_qm(atom_1, atom_2, 0.96, 0.02), zip(atoms_1, atoms_2))
        #print "BL DEBUG:: all true", all_true
        # A bit of a windy one...
        msg = ""
        if all_false:
            msg = "  Oops! bond length not within tolerance: %s\n  Hydrogen names mangled from PDB %s %s" %(all_false[0][0],
                                                                                                            all_false[0][1], all_false[0][2])
        self.failUnless(len(all_false) == 0, msg)


    def test32_1(self):
        """Correct matching dictionary names from test name"""

        # new dictionary
        ls = matching_compound_names_from_dictionary("gua", 0)
        # too fragile
        #self.failUnless(len(ls) == 153 or \  # new dictionary
        #                len(ls) == 63)       # old dicionary
        self.failUnless(len(ls) > 60,
                        "   found %s matching names: %s" %(len(ls), ls))


    def test33_0(self):
        """Update monomer restraints"""

        atom_pair = [" CB ", " CG "]
        m = monomer_restraints("TYR")
        self.failUnless(m, "   update bond restraints - no momomer restraints")

        # let's delete both ways
        #
        n_t = strip_bond_from_restraints(atom_pair, m)
        atom_pair.reverse()
        n = strip_bond_from_restraints(atom_pair, n_t)
        set_monomer_restraints("TYR", n)

        imol = new_molecule_by_atom_selection(imol_rnase, "//A/30")

        atom_1 = get_atom(imol, "A", 30, "", " CB ")
        atom_2 = get_atom(imol, "A", 30, "", " CG ")

        with_auto_accept([refine_zone, imol, "A", 30, 30, ""])

        atom_1 = get_atom(imol, "A", 30, "", " CB ")
        atom_2 = get_atom(imol, "A", 30, "", " CG ")
        print "  Bond-length: ", bond_length(atom_1[2], atom_2[2])

        self.failUnless(bond_length_within_tolerance_qm(atom_1, atom_2, 2.8, 0.6),
                        "fail 2.8 tolerance test")
        print "   pass intermediate 2.8 tolerance test"
        set_monomer_restraints("TYR", m)

        with_auto_accept([refine_zone, imol, "A", 30, 30, ""])

        atom_1 = get_atom(imol, "A", 30, "", " CB ")
        atom_2 = get_atom(imol, "A", 30, "", " CG ")

        post_set_plane_restraints = monomer_restraints("TYR")["_chem_comp_plane_atom"]

        atom = post_set_plane_restraints[0][1][0]
        self.failUnless(atom == " CB ", "FAIL plane atom %s" %atom)
        print "  OK plane atom ", atom

        print "  Bond-length: ", bond_length(atom_1[2], atom_2[2])
        self.failUnless(bond_length_within_tolerance_qm(atom_1, atom_2, 1.512, 0.04),
                        "fail 1.512 tolerance test")


    def test33_1(self):
        """Write mmCIF restraints correctly"""

        # not a very good test, but better than nothing and it would
        # fail on cif writing prior to today 20100223?!

        if os.path.isfile("coot-test-ala.cif"):
            os.remove("coot-test-ala.cif")

        write_restraints_cif_dictionary("ALA", "coot-test-ala.cif")
        self.failUnless(os.path.isfile("coot-test-ala.cif"))
        fin = open("coot-test-ala.cif", 'r')
        lines = fin.readlines()
        fin.close()
        n_found = 0
        for line in lines:
            if "data_comp_list" in line:
                n_found = 1
            if "data_comp_ALA" in line:
                n_found += 1
        self.failUnlessEqual(n_found, 2)
        
    def test34_0(self):
        """Refinement OK with zero bond esd"""

        # return the restraints as passed (in r_list) except for a bond
        # restraint atom_1 to atom_2 set the esd to 0.0 (new_dist is
        # passed but not useful).
        def zero_bond_restraint(r_list, atom_1, atom_2, new_dist):

            bond_list = r_list["_chem_comp_bond"]
            for restraint in bond_list:
                if (restraint[0] == atom_1 and
                    restraint[1] == atom_2):
                    bond_list[bond_list.index(restraint)] = [atom_1,
                                                             atom_2,
                                                             new_dist,
                                                             0.0]
            r_list["_chem_comp_bond"] = bond_list
            return r_list

        # main line
        #
        import os
        from types import DictType
        imol = unittest_pdb("monomer-ACT.pdb")
        read_cif_dictionary(os.path.join(unittest_data_dir,
                                         "libcheck_ACT.cif"))
        self.failUnless(valid_model_molecule_qm(imol),
                        "   bad molecule from ACT from greg data dir")

        r = monomer_restraints("ACT")
        self.failUnless(type(r) is DictType,
                        "   ACT restraints are False")

        r_2 = zero_bond_restraint(r, " O  ", " C  ", 1.25)
        #print "r:  ", r
        #print "r2: ", r_2
        self.failIf(r_2 == {}, "   null modified restraints")

        mr_set = set_monomer_restraints("ACT", r_2)
        self.failUnless(mr_set, "   set restraints fail")

        r_3 = monomer_restraints("ACT")
        self.failUnless(r_3, "   get modified restraints fail")

        with_auto_accept([refine_zone, imol, "A", 1, 1, ""])

        # didn't crash?!

# enable when we are in 07
# simple-minded Hydrogen-based filtered (not topological).
#    def test34_1(self):
#        """filter out multi-hydrogen chiral centres"""
#
#        global unittest_data_dir
#        import os
#        imol = unittest_pdb("chiral-test-1.pdb")
#        self.failUnless(valid_model_molecule_qm(imol),
#                        "chiral-test-1.pdb not found")
#        read_cif_dictionary(os.path.join(unittest_data_dir,
#                                         "chiral-test-1-dict.cif"))
#        r = monomer_restraints("DRG")
#        chirals = r["_chem_comp_chir"]
#        self.failIf(len(chirals) == 2)  # ? FIXME

    def test35_0(self):
        """Change Chain IDs and Chain Sorting"""

        def chains_in_order_qm(chain_list):
            ref_chain = ""
            for chain in chain_list:
                if (chain < ref_chain):
                    print "ERROR:: %s was less than %s in %s" %(chain, ref_chain, chain_list)
                    return False
                else:
                    ref_chain = chain
            return True

        imol = unittest_pdb("tutorial-modern.pdb")

        change_chain_id(imol, "A", "D", 0,  0,  0)
        change_chain_id(imol, "B", "N", 1, 30, 38)
        change_chain_id(imol, "B", "L", 1, 40, 49)
        change_chain_id(imol, "B", "J", 1, 50, 59)
        change_chain_id(imol, "B", "Z", 1, 20, 28)
        change_chain_id(imol, "B", "G", 1, 60, 70)
        change_chain_id(imol, "B", "F", 1, 70, 80)
        change_chain_id(imol, "B", "E", 1, 80, 90)

        sort_chains(imol)

        c = chain_ids(imol)
        #print "BL DEBUG:: chains_in_order_qm", chains_in_order_qm(c)
        self.failUnless(chains_in_order_qm(c))


    def test35_1(self):
        """Chain-ids in links change also on change chain id"""

        imol = unittest_pdb("tutorial-modern.pdb")

        spec_1 = ["B", 42, "", " N  ", ""]
        spec_2 = ["B", 62, "", " O  ", ""]

        make_link(imol, spec_1, spec_2, "test-link", 2.2)

        li_1 = link_info(imol)

        change_chain_id(imol, "B", "C", 0, 0, 0)

        li_2 = link_info(imol)

        # li-1 should not contain C and should contain B
        # li-2 should not contain B and should contain C

        ch_B_1 = li_1[0][1][1] # before
        ch_B_2 = li_1[0][2][1] # 
        ch_A_1 = li_2[0][1][1] # after
        ch_A_2 = li_2[0][2][1] # 

        self.failUnless(all([ch_B_1 == "B",
                             ch_B_2 == "B",
                             ch_A_1 == "C",
                             ch_A_2 == "C"]))

        
    def test36_0(self):
        """Replace Fragment"""

        atom_sel_str = "//A70-80"
        imol_rnase_copy = copy_molecule(imol_rnase)
        imol = new_molecule_by_atom_selection(imol_rnase, atom_sel_str)
        self.failUnless(valid_model_molecule_qm(imol))

        translate_molecule_by(imol, 11, 12, 13)
        reference_res = residue_info(imol_rnase, "A", 75, "")
        moved_res     = residue_info(imol,       "A", 75, "")
        test_res     = residue_info(imol,       "A", 74, "")

        replace_fragment(imol_rnase_copy, imol, atom_sel_str)
        replaced_res = residue_info(imol_rnase_copy, "A", 75, "")

        # now replace-res should match reference-res.
        # the atoms of moved-res should be 20+ A away from both.

        self.failUnless(all(map(atoms_match_qm, moved_res, replaced_res)),
                        "   moved-res and replaced-res do not match")

        self.failIf(all(map(atoms_match_qm, moved_res, reference_res)),
                    "   fail - moved-res and replaced-res Match!")

        print "   distances: ", map(atom_distance, reference_res, replaced_res)
        self.failUnless(all(map(lambda d: d > 20, map(atom_distance, reference_res, replaced_res))))



    def test37_0(self):
        """Residues in Region of Residue"""

        for dist, n_neighbours in zip([4,0], [6,0]):
            rs = residues_near_residue(imol_rnase, ["A", 40, ""], dist)
            self.failUnless(len(rs) == n_neighbours,
                            "wrong number of neighbours %s %s" %(len(rs), rs))
            print "   found %s neighbours %s" %(len(rs), rs)


    def test38_0(self):
        """Residues in region of a point"""

        imol = unittest_pdb("tutorial-modern.pdb")
        pt = [43.838, 0.734, 13.811] # CA 47 A
        residues = residues_near_position(imol, pt, 2)
        self.failUnless(len(residues) == 1,
                        "  Fail, got residues: %s" %residues)
        self.failUnless(residues[0] == [True, "A", 47, ""],
                        "  Fail 2, got residues: %s" %residues)
        residues_2 = residues_near_position(imol, pt, 4)
        self.failUnless(len(residues_2) == 3,
                        "  Fail 3, got residues-2: %s" %residues_2) #its neighbours too.


    def test39_0(self):
        """Empty molecule on type selection"""

        global imol_rnase
        imol1 = new_molecule_by_residue_type_selection(imol_rnase, "TRP")
        self.failUnlessEqual(imol1, -1,
                             "failed on empty selection 1 gives not imol -1")
        imol2 = new_molecule_by_residue_type_selection(imol_rnase, "TRP")
        self.failUnlessEqual(imol2, -1,
                             "failed on empty selection 2 gives not imol -1")


    def test40_0(self):
        """Set Rotamer"""

        chain_id = "A"
        resno = 51

        n_rot = n_rotamers(-1, "ZZ", 45, "")
        self.failUnless(n_rot == -1)

        n_rot = n_rotamers(imol_rnase, "Z", 45, "")
        self.failUnless(n_rot == -1)

        residue_pre = residue_info(imol_rnase, chain_id, resno, "")

        # note that the rotamer number is 0-indexed (unlike the rotamer
        # number dialog)
        set_residue_to_rotamer_number(imol_rnase, chain_id, resno, "", "", 1)

        residue_post = residue_info(imol_rnase, chain_id, resno, "")

        self.failUnless(len(residue_pre) == len(residue_post))

        # average dist should be > 0.1 and < 0.3.
        #
        dists = map(atom_distance, residue_pre, residue_post)
        #print "BL DEBUG:: dists", dists
        self.failUnless(all(map(lambda d: d <= 0.6, dists)))
        self.failUnless(all(map(lambda d: d >= 0.0, dists)))


    def test41_0(self):
        """Rotamer names and scores are correct"""

        imol = unittest_pdb("tutorial-modern.pdb")
        turn_off_backup(imol)  # easier for long function that using with_no_backup
        residue_attributes = ["A", 28, ""]
        for rotamer_number, correct_name, correct_prob \
            in zip(range(5),
                   # the spaces at the end of the name are in the Complete_rotamer_lib.csv
                   # and hence in richardson-rotamers.cc.  C'est la vie.
                   ["m-85", "t80", "p90", "m -30 ", "m -30 "],
                   [100, 90.16684, 50.707787, 21.423154, 21.423154]):
            set_residue_to_rotamer_number(imol, *(residue_attributes + ["", rotamer_number]))
            rotamer_name = get_rotamer_name(imol, *residue_attributes)
            rotamer_prob = rotamer_score(imol, *(residue_attributes + [""]))
            print "   Rotamer %s : %s %s" %(rotamer_number, rotamer_name, rotamer_prob)
            self.failUnlessAlmostEqual(rotamer_prob, correct_prob, 3,
                                       "fail on rotamer probability: result: %s %s" %(rotamer_prob, correct_prob)),
            self.failUnless(rotamer_name == correct_name)

        turn_on_backup(imol)


    def test42_0(self):
        """Align and mutate a model with deletions"""

        def residue_in_molecule_qm(imol, chain_id, resno, ins_code):
            r = residue_info(imol, chain_id, resno, ins_code)
            if (r):
                return True
            else:
                return False

        # in this PDB file 60 and 61 have been deleted. Relative to the
        # where we want to be (tutorial-modern.pdb, say) 62 to 93 have
        # been moved to 60 to 91
        #
        imol = unittest_pdb("rnase-A-needs-an-insertion.pdb")
        if self.skip_test(not os.path.isfile(rnase_seq),
                          "file %s not found, skipping test" %rnase_seq):
            return
        rnase_seq_string = file2string(rnase_seq)
        self.failUnless(valid_model_molecule_qm(imol),
                        "   Missing file rnase-A-needs-an-insertion.pdb")
        renumber = 1

        set_alignment_gap_and_space_penalty(-3.0, -0.5)
        align_and_mutate(imol, "A", rnase_seq_string, renumber)
        write_pdb_file(imol, "mutated.pdb")

        ls = [[False, imol, "A",  1, ""],
              [True,  imol, "A",  4, ""],
              [True,  imol, "A", 59, ""],
              [False, imol, "A", 60, ""],
              [False, imol, "A", 61, ""],
              [True,  imol, "A", 93, ""],
              [False, imol, "A", 94, ""]]

        def compare_res_spec(res_info):
            residue_spec = res_info[1:len(res_info)]
            expected_status = res_info[0]
            print "    :::::", residue_spec,
            residue_in_molecule_qm(*residue_spec), expected_status
            if (residue_in_molecule_qm(*residue_spec) == expected_status):
                return True
            else:
                return False

        results = map(lambda res: compare_res_spec(res), ls)
        print "results", results
        self.failUnless(all(results))


    def test42_1(self):
        """Renumber residues should be in seqnum order"""

        imol = unittest_pdb("tutorial-modern.pdb")
        self.failUnless(valid_model_molecule_qm(imol),
                        " ERROR:: tutorial-modern.pdb not found")
        renumber_residue_range(imol, "A", 10, 20, 100)
        iser = seqnum_from_serial_number(imol, "A", 90)
        self.failUnlessEqual(iser, 118)

        # and again this time adding residue to the front of the molecule
        imol = unittest_pdb("tutorial-modern.pdb")
        renumber_residue_range(imol, "A", 10, 20, -100)
        iser = seqnum_from_serial_number(imol, "A", 0)
        self.failUnlessEqual(iser, -90)
        
        
    def test43_0(self):
        """Autofit Rotamer on Residue with Insertion codes"""

        # we need to check that H52 LEU and H53 GLY do not move and H52A does move

        def centre_atoms(mat):
            tm = transpose_mat(map(lambda x: x[2], mat))
            centre = map(lambda ls: sum(ls)/len(ls), tm)
            return centre

        #
        imol = unittest_pdb("pdb3hfl.ent")
        mtz_file_name = os.path.join(unittest_data_dir, "3hfl_sigmaa.mtz")
        imol_map = make_and_draw_map(mtz_file_name,
                                     "2FOFCWT", "PH2FOFCWT", "", 0, 0)

        self.failUnless(valid_model_molecule_qm(imol),
                        " ERROR:: pdb3hfl.ent not found")

        self.failUnless(valid_map_molecule_qm(imol_map),
                        "   no map from 3hfl_sigmaa.mtz")

        leu_atoms_1 = residue_info(imol, "H", 52, "")
        leu_resname = residue_name(imol, "H", 52, "")
        gly_atoms_1 = residue_info(imol, "H", 53, "")
        gly_resname = residue_name(imol, "H", 53, "")
        pro_atoms_1 = residue_info(imol, "H", 52, "A")
        pro_resname = residue_name(imol, "H", 52, "A")

        # First check that the residue names are correct
        self.failUnless((leu_resname == "LEU" and
                         gly_resname == "GLY" and
                         pro_resname == "PRO"),
                        "  failure of residues names: %s %s %s" \
                        %(leu_resname, gly_resname, pro_resname))

        # OK, what are the centre points of these residues?
        leu_centre_1 = centre_atoms(leu_atoms_1)
        gly_centre_1 = centre_atoms(gly_atoms_1)
        pro_centre_1 = centre_atoms(pro_atoms_1)

        auto_fit_best_rotamer(52, "", "A", "H", imol, imol_map, 0, 0)

        # OK, what are the centre points of these residues?
        leu_atoms_2  = residue_info(imol, "H", 52, "")
        gly_atoms_2  = residue_info(imol, "H", 53, "")
        pro_atoms_2  = residue_info(imol, "H", 52, "A")
        leu_centre_2 = centre_atoms(leu_atoms_2)
        gly_centre_2 = centre_atoms(gly_atoms_2)
        pro_centre_2 = centre_atoms(pro_atoms_2)

        d_leu = pos_diff(leu_centre_1, leu_centre_2)
        d_gly = pos_diff(gly_centre_1, gly_centre_2)
        d_pro = pos_diff(pro_centre_1, pro_centre_2)

        self.failUnlessAlmostEqual(d_leu, 0.0, 1,
                                   "   Failure: LEU 52 moved")

        self.failUnlessAlmostEqual(d_gly, 0.0, 1,
                                   "   Failure: GLY 53 moved")

        self.failIf(d_pro < 0.05,
                    """   Failure: PRO 52A not moved enough %s\n 
                    PRO-atoms 1: %s\n
                    PRO-atoms 2: %s\n""" %(d_pro, pro_atoms_1, pro_atoms_2,))
        # in theory (greg-tests) there should extra output for each atom
        # in pro-x

        # rotamer 4 is out of range.
        set_residue_to_rotamer_number(imol, "H", 52, "A", "", 4) # crash!!?


    # new version of tests
    def test44_0(self):
        """RNA base has correct residue type after mutation"""

        def test_vs(rna_mol, base_name):

            previous_name = residue_name(rna_mol, "A", 2, "")
            success = mutate_base(rna_mol, "A", 2, "", base_name)
            self.failUnless(success == 1, "  mutation fail!")
            rn = residue_name(rna_mol, "A", 2, "")
            print "  mutated base to type %s - was %s" %(rn, previous_name)
            return (rn == base_name)

        # main line
        #   
        rna_mol = ideal_nucleic_acid("RNA", "A", 0, "GACUCUAG")
        res_1 = test_vs(rna_mol, "C")
        self.failUnless(res_1, "  incorrect base! (default names)")

        set_convert_to_v2_atom_names(1)

        rna_mol_old_names = ideal_nucleic_acid("RNA", "A", 0, "GACUCUAG")
        res_2 = test_vs(rna_mol_old_names, "Cr")

        # back to normal
        set_convert_to_v2_atom_names(0)
        
        self.failUnless(res_2, "  incorrect base! (old names)")


    def test45_0(self):
        """DNA bases are the correct residue type after mutation"""

        def correct_base_type_qm(rna_mol, target_base_type):
            success = mutate_base(rna_mol, "A", 2, "", target_base_type)
            self.failUnless(success == 1, "  DNA base mutation fail!")
            rn = residue_name(rna_mol, "A", 2, "")
            print "  mutated base to type", rn
            return (rn == target_base_type)

        # main line
        rna_mol = ideal_nucleic_acid("DNA", "A", 0, "GACTCTAG")

        res = map(lambda base: correct_base_type_qm(rna_mol, base),
                  ["DC", "DG", "DA", "DT"])
        self.failUnless(all(res))

        set_convert_to_v2_atom_names(1)
        rna_mol = ideal_nucleic_acid("DNA", "A", 0, "GACTCTAG")
        
        res = map(lambda base: correct_base_type_qm(rna_mol, base),
                  ["Cd", "Gd", "Ad", "Td"])
        set_convert_to_v2_atom_names(0)
        self.failUnless(all(res))


    def test46_0(self):
        """SegIDs are correct after mutate"""

        # main line
        global imol_rnase
        imol = copy_molecule(imol_rnase)

        self.failUnless(valid_model_molecule_qm(imol))

        atoms = residue_info(imol, "A", 32, "")

        self.failUnless(atoms_have_correct_seg_id_qm(atoms, ""),
                        "wrong seg-id %s should be %s" %(atoms, ""))

        # now convert that residue to segid "A"
        #
        attribs = map(lambda atom: [imol, "A", 32, "",
                                    atom[0][0],
                                    atom[0][1],
                                    "segid", "A"],
                      residue_info(imol, "A", 32, ""))
        set_atom_attributes(attribs)

        atoms = residue_info(imol, "A", 32, "")
        self.failUnless(atoms_have_correct_seg_id_qm(atoms, "A"),
                        "wrong seg-id %s should be %s" %(atoms, "A"))

        # now let's do the mutation
        mutate(imol, "A", 32, "", "LYS")

        rn = residue_name(imol, "A", 32, "")
        self.failUnless(rn == "LYS",
                        "  Wrong residue name after mutate %s" %rn)

        atoms = residue_info(imol, "A", 32, "")
        self.failUnless(atoms_have_correct_seg_id_qm(atoms, "A"),
                        "wrong seg-id %s should be %s" %(atoms, "A"))


    def test47_0(self):
        """TER on water chain is removed on adding a water by hand"""

        imol = unittest_pdb("some-waters-with-ter.pdb")

        self.failUnless(valid_model_molecule_qm(imol),
                        "bad read of some-waters-with-ter.pdb")

        set_rotation_centre(3, 4, 5)
        set_pointer_atom_molecule(imol)
        place_typed_atom_at_pointer("Water")

        # OK let's write out that molecule now as a PDB file and look
        # to see if the PDB file contains a TER record - it should
        # not.
        #
        opdb = "tmp-with-new-water.pdb"
        write_pdb_file(imol, opdb)
        fin = open(opdb, 'r')
        lines = fin.readlines()
        fin.close()
        found_ter = False
        for line in lines:
            self.failIf("TER" in line,
                        "   TER card found: %s" %line)


    def test48_0(self):
        """TER on water chain is removed on adding waters automatically"""

        global imol_rnase_map
        imol_model = unittest_pdb("tm+some-waters.pdb")
        self.failUnless(valid_model_molecule_qm(imol_model),
                        "tm+some-waters.pdb not found")
        find_waters(imol_rnase_map, imol_model, 0, 0.2, 0)

        # OK let's write out that molecule now as a PDB file and look
        # to see if the PDB file contains a TER record - it should
        # not.
        #
        opdb = "auto-waters.pdb"
        write_pdb_file(imol_model, opdb)
        fin = open(opdb, 'r')
        lines = fin.readlines()
        fin.close()
        found_ter = False
        for line in lines:
            self.failIf("TER" in line,
                        "   TER card found: %s" %line)


    def test48_1(self):
        """Adding atoms to Many-Chained Molecule"""

        imol = read_pdb(rnase_pdb())
        set_pointer_atom_molecule(imol)
        for i in range(100):
            place_typed_atom_at_pointer("Mg")
        # doesnt crash...
        
            
    def test49_0(self):
        """Arrange waters round protein"""

        imol = unittest_pdb("water-test-no-cell.pdb")
        self.failUnless(valid_model_molecule_qm(imol),
                        "ERROR:: water-test-no-cell not found")

        status = move_waters_to_around_protein(imol)
        self.failUnless(status == 0,
                        "ERROR:: failure with water-test-no-cell")

        imol = unittest_pdb("pathological-water-test.pdb")
        self.failUnless(valid_model_molecule_qm(imol),
                        "ERROR:: pathological-waters pdb not found")

        status = move_waters_to_around_protein(imol)
        write_pdb_file(imol, "waters-moved-failure.pdb")
        self.failUnless(status == 181,
                        "ERROR:: failure with pathological-waters moved %s" %(status))

        v = max_water_distance(imol)

        self.failIf(v > 5.0,
                    "ERROR:: failure to move waters close %s" %v)


    def test50_0(self):
        """Correct Segid After Add Terminal Residue"""

        global imol_rnase
        imol = copy_molecule(imol_rnase)
        imol_map = make_and_draw_map(rnase_mtz(), "FWT", "PHWT", "", 0, 0)

        # now convert that residue to segid "A"
        #
        attribs = map(lambda atom: [imol, "A", 93, "",
                                    atom[0][0],
                                    atom[0][1],
                                    "segid", "A"],
                      residue_info(imol, "A", 93, ""))
        set_atom_attributes(attribs)

        add_terminal_residue(imol, "A", 93, "ALA", 1)

        atoms = residue_info(imol, "A", 94, "")
        self.failUnless(atoms_have_correct_seg_id_qm(atoms, "A"),
                        "wrong seg-id %s should be %s" %(atoms, "A"))


    def test51_0(self):
        """Correct Segid after NCS residue range copy"""

        def convert_residue_seg_ids(imol, chain_id, resno, seg_id):
            attribs = map(lambda atom: [imol, chain_id, resno, "",
                                        atom[0][0],
                                        atom[0][1],
                                        "segid", seg_id],
                          residue_info(imol, chain_id, resno, ""))
            set_atom_attributes(attribs)

        # main line
        global imol_rnase
        imol = copy_molecule(imol_rnase)

        # convert a residue range
        for resno in range(20, 30):
            convert_residue_seg_ids(imol, "A", resno, "X")

        # NCS copy
        make_ncs_ghosts_maybe(imol)
        copy_residue_range_from_ncs_master_to_others(imol, "A", 20, 30)

        # does B have segid X for residues 20-30?
        for resno in range(20, 30):
            atoms = residue_info(imol, "A", resno, "")
            self.failUnless(atoms_have_correct_seg_id_qm(atoms, "X"),
                            "wrong seg-id %s should be %s" %(atoms, "X"))



    def test52_0(self):
        """Merge Water Chains"""

        def create_water_chain(imol, from_chain_id, chain_id, n_waters,
                               offset, prev_offset):
            for n in range(n_waters):
                place_typed_atom_at_pointer("Water")
                # move the centre of the screen
                rc = rotation_centre()
                set_rotation_centre(rc[0] + 2., rc[1], rc[2])
            change_chain_id(imol, from_chain_id, chain_id, 1,
                            prev_offset + 1, n_waters + prev_offset)
            renumber_residue_range(imol, chain_id, prev_offset + 1,
                                   prev_offset + 5, offset)

        # main line
        imol = unittest_pdb("tutorial-modern.pdb")
        # with_no_backup....
        turn_off_backup(imol)
        set_pointer_atom_molecule(imol)

        # first create a few chains of waters
        create_water_chain(imol, "C", "D", 5, 10,  0)
        create_water_chain(imol, "D", "E", 5,  5, 10)
        create_water_chain(imol, "D", "F", 5,  2, 15)

        # OK, we have set up a molecule, now let's test it:
        #
        merge_solvent_chains(imol)
        turn_on_backup(imol)

        # Test the result:
        #
        nc = n_chains(imol)
        self.failUnless(nc == 3,
                        "  wrong number of chains %s" %nc)
        # There should be 15 waters in the last chain
        solvent_chain_id = chain_id(imol, 2)
        n_res = chain_n_residues(solvent_chain_id, imol)
        self.failUnless(n_res == 15,
                        "  wrong number of residues %s" %n_res)
        r1  = seqnum_from_serial_number(imol, "D", 0)
        r15 = seqnum_from_serial_number(imol, "D", 14)
        self.failUnless(r1  == 1,
                        "  wrong residue number r1 %s"  %r1)
        self.failUnless(r15 == 15,
                        "  wrong residue number r15 %s" %r15)


    def test52_1(self):
        """Consolidated merge"""

        imol = unittest_pdb("pdb1hvv.ent")
        imol_lig_1 = unittest_pdb("monomer-ACT.pdb")
        imol_lig_2 = unittest_pdb("monomer-NPO.pdb")

        print "-------- starting chain list -----------"
        print chain_ids(imol)

        merge_molecules([imol_lig_1, imol_lig_2], imol)

        imol_symm_copy = new_molecule_by_symop(imol, "-X,-X+Y,-Z+1/3", 0, 0, 0)

        self.failUnless(valid_model_molecule_qm(imol_symm_copy),
                        "Symm molecule problem")

        merge_molecules([imol_symm_copy], imol)

        chain_list = chain_ids(imol)

        print chain_list

        write_pdb_file(imol, "sym-merged.pdb")  # why?

        self.failUnlessEqual(chain_list,
                             ["A", "B", "C", "D", "",  # original
                              "E", "F",                # merged ligs
                              "G", "H", "I", "J", "K"],  # protein chain copies
                             "List did not match")
        

    def test52_2(self):
        """Test for good chain ids after a merge"""

        imol = unittest_pdb("tutorial-modern.pdb")

        change_chain_id(imol, "A", "AAA", 0, 0, 0)

        imol_new = new_molecule_by_atom_selection(imol, "//B/1-90")

        change_chain_id(imol_new, "B", "B-chain", 0, 0, 0)

        merge_molecules([imol_new], imol)

        chids = chain_ids(imol)

        print "chain_ids", chids

        # should be ["AAA", "B", "B-chain"]

        self.failUnless(len(chids) == 3)

        self.failUnless(chids[2] == "B-chain")

        # now Wolfram Tempel test: multi-char chain matcher needs
        # prefix, not new single letter

        change_chain_id(imol_new, "B-chain", "AAA", 0, 0, 0)

        merge_molecules([imol_new], imol)

        chids_2 = chain_ids(imol)

        print "--- chain-ids:", chids_2

        self.failUnless(len(chids_2) == 4)

        self.failUnless(chids_2[3] == "AAA2")
        

    def test53_0(self):
        """LSQ by atom"""

        global imol_rnase

        def make_spec_ref(atom_name):
            return ["A", 35, "", atom_name, ""]

        def make_spec_mov(atom_name):
            return ["B", 35, "", atom_name, ""]

        clear_lsq_matches()
        imol_1= copy_molecule(imol_rnase)
        imol_2= copy_molecule(imol_rnase)

        spec_refs = map(make_spec_ref, [" CG2", " CG1", " CB ", " CA "])
        spec_movs = map(make_spec_mov, [" CG2", " CG1", " CB ", " CA "])
        map(add_lsq_atom_pair, spec_refs, spec_movs)

        result = apply_lsq_matches(imol_1, imol_2)
        self.failUnless(result, "Bad match")

        c_1 = get_atom(imol_1, "A", 35, "", " C  ")
        c_2 = get_atom(imol_2, "B", 35, "", " C  ")
        b = bond_length_from_atoms(c_1, c_2)

        self.failUnless(bond_length_within_tolerance_qm(c_1, c_2, 0.0, 0.2))


    def test53_1(self):
        """LSQing changes the space-group and cell to that of the reference molecule"""

        imol_mov = unittest_pdb("tutorial-modern.pdb")
        imol_ref = unittest_pdb("pdb1py3.ent")

        sg_mov_orig = show_spacegroup(imol_mov)
        sg_ref_orig = show_spacegroup(imol_ref)
        cell_mov_orig = cell(imol_mov)
        cell_ref_orig = cell(imol_ref)

        clear_lsq_matches()
        add_lsq_match(10, 50, "A", 8, 48, "B", 1)
        rtop = apply_lsq_matches(imol_ref, imol_mov)
        sg_mov_curr = show_spacegroup(imol_mov)
        cell_mov_curr = cell(imol_mov)

        self.failUnless(sg_mov_curr == sg_ref_orig,
                        "   fail on matching spacegroups: %s and %s" \
                        %(sg_mov_curr, sg_ref_orig))

        self.failUnless(all(map(close_float_qm, cell_ref_orig, cell_mov_curr)),
                        "   fail on matching cells: %s and %s" \
                        %(cell_ref_orig, cell_mov_curr))


    def test53_2(self):
        """set_residue_name sets the correct residue"""

        imol = unittest_pdb("tutorial-modern.pdb")
        set_residue_name(imol, "A", 37, "", "FRE")

        # There should be only one residue with that residue type and it
        # should be A37.
        #
        specs = fit_protein_make_specs(imol, 'all-chains')
        residue_with_name = get_residues_in_molecule_of_type(imol, "FRE")

        self.failUnless(len(residue_with_name) == 1)
        self.failUnless(spec_match_qm(residue_with_name[0], [imol, "A", 37, ""]))


    def test53_3(self):
        """fit_protein_make_specs makes all specs"""

        imol = unittest_pdb("tutorial-modern.pdb")
        specs = fit_protein_make_specs(imol, 'all-chains')
        print "   specs:", len(specs), specs
        self.failUnless(len(specs) == 189)
   
        

    def test54_0(self):
        """Phosphate distance in pucker analysis is sane"""

        imol = unittest_pdb("2goz-manip.pdb")
        self.failUnless(valid_model_molecule_qm(imol),
                        "failed to find 2goz-manip.pdb")

        pi = pucker_info(imol, ["B", 14, ""], 0)
        phosphate_distance = pi[0]
        self.failIf(phosphate_distance > 0.2,
                    "  Bad phosphate distance on 14 %s" %phosphate_distance)

        pi = pucker_info(imol, ["B", 15, ""], 0)
        phosphate_distance = pi[0]
        self.failIf(phosphate_distance < 2.0,
                    "  Bad phosphate distance on 15 %s" %phosphate_distance)


    # replacement of Ramachandran refinement and just testing
    # post_manipulation_hook_py (out of context)
    def test60_0(self):
        """post_manipulation_hook_py test (replaces Hundreds of Ramachandran refinements)"""

        # doesnt test for anything just crash
        n_post = 10000
        for i in range(n_post):
            safe_python_command_with_return("post_manipulation_script")


    def test61_0(self):
        """Read/write gz coordinate files"""

        # this is mainly an mmdb test for windows

        # first write a gz file
        gz_state = write_pdb_file(imol_rnase, "rnase_zip_test.pdb.gz")
        self.failIf(gz_state == 1)
        self.failUnless(os.path.isfile("rnase_zip_test.pdb.gz"))

        # now unzip and read the file
        gz_imol = handle_read_draw_molecule("rnase_zip_test.pdb.gz")
        self.failUnless(valid_model_molecule_qm(gz_imol))

        
    def test62_0(self):
        """Fix for Oliver Clarke fit by atom selection bug"""

        imol_rnase = unittest_pdb("tutorial-modern.pdb")
        imol_map = make_and_draw_map(rnase_mtz(), "FWT", "PHWT", "", 0, 0)

        set_imol_refinement_map(imol_map)

        with AutoAccept():
            rigid_body_refine_by_atom_selection(imol_rnase, "//B")

        # did we get silly superposition?
        # test that the A-chain atom is close to where it should be

        atom = get_atom(imol_rnase, "A", 43, "", " CA ", "")
        self.failUnless(isinstance(atom, list), "Failure to extract atom")
        atom_pos = atom[2]
        bl = bond_length(atom_pos, [46.4, 11.6, 12.1])
        print "bl:", bl
        self.failIf(bl > 1.0, "Fail: moved atom %s" %bl)


    def test999_0(self):
        """Renumber residue range without overwriting coordinates."""

        # This tests that we
        imol = unittest_pdb("tutorial-modern.pdb")
        self.failUnless(valid_model_molecule_qm(imol))
        self.failUnless(renumber_residue_range(imol, "A", 10, 20, -55) == 1)
        self.failUnless(renumber_residue_range(imol, "A", 90, 93, 10) == 1)
        self.failUnless(renumber_residue_range(imol, "A", 89, 91, 1) == 1)
        self.failUnless(renumber_residue_range(imol, "A", 80, 91, 12) == 0)
        self.failUnless(renumber_residue_range(imol, "A", -100, 200, 9) == 1)
        
        
