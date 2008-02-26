# Copyright 2007, 2008 by The University of York
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

rnase_pdb = os.path.join(unittest_data_dir, "tutorial-modern.pdb")
rnase_mtz = os.path.join(unittest_data_dir, "rnasa-1.8-all_refmac1.mtz")
terminal_residue_test_pdb = os.path.join(unittest_data_dir, "tutorial-add-terminal-0-test.pdb")
base_imol = graphics_n_molecules()

global have_ccp4_qm, imol_rnase
have_ccp4_qm = False
imol_rnase = -1
imol_ligand = -1
imol_terminal_residue_test = -1

horne_works_cif = os.path.join(unittest_data_dir, "lib-B3A.cif")
horne_cif = os.path.join(unittest_data_dir, "lib-both.cif")
horne_pdb = os.path.join(unittest_data_dir, "coords-B3A.pdb")

# CCP4 is set up? If so, set have-ccp4? True

try:
	#global have_ccp4_qm
	ccp4_master = os.getenv("CCP4_MASTER")
	if (os.path.isfile(ccp4_master)):
		have_ccp4_qm = True
except:
	print "BL INFO:: Dont have CCP4 master"

class PdbMtzTestFunctions(unittest.TestCase):

    # tests are executed alphanumerical, so we shall give them number, rather than
    # names. We add a 0 in the end to give space for later addition
    def test01_0(self):
        """Close bad molecule"""
        close_molecule(-2)

    def test02_0(self):
	"""Read coordinates test"""
	imol = read_pdb(rnase_pdb)
	global imol_rnase
	imol_rnase = imol
	self.failUnless(valid_model_molecule_qm(imol))

    def test03_0(self):
	    """New molecule from bogus molecule"""
	    pre_n_molecules = graphics_n_molecules()
	    new_mol = new_molecule_by_atom_selection(-5, "//A/0")
	    post_n_molecules = graphics_n_molecules()
	    print "BL DEBUG:: new_mol, pre and post", new_mol, pre_n_molecules, post_n_molecules
	    self.failUnlessEqual(new_mol, -1 ,"fail on non-matching n molecules (-1)")
	    self.failUnlessEqual(pre_n_molecules, post_n_molecules, "fail on non-matching n molecules (=)")

    def test04_0(self):
	    """New molecule from bogus atom selection"""
	    global imol_rnase
	    pre_n_molecules = graphics_n_molecules()
	    new_molecule = new_molecule_by_atom_selection(imol_rnase, "//A/0")
	    post_n_molecules = graphics_n_molecules()
	    # what should happen if there are no atoms in the new-mol?
	    print "   INFO:: pre_n_molecules %s   post_n_molecules %s" %(pre_n_molecules, post_n_molecules)
	    self.failUnlessEqual(new_molecule, -1 ,"fail on non-matching n molecules (-1)")
	    self.failUnlessEqual(pre_n_molecules, post_n_molecules, "fail on non-matching n molecules (=)")

    def test05_0(self):
	    """Read a bogus map"""
	    pre_n_molecules = graphics_n_molecules()
	    imol = handle_read_ccp4_map("bogus.map", 0)
	    print "BL DEBUG:: bogus map imol", imol
	    self.failUnlessEqual(imol, -1, "bogus ccp4 map returns wrong molecule number")
	    now_n_molecules = graphics_n_molecules()
	    self.failUnlessEqual(now_n_molecules, pre_n_molecules, "bogus ccp4 map creates extra map %s %s " %(pre_n_molecules, now_n_molecules))
	    
    def test06_0(self):
	    """Read MTZ test"""
	    
	    # bogus map test
	    pre_n_molecules = graphics_n_molecules()
	    imol_map = make_and_draw_map("bogus.mtz", "FWT", "PHWT", "", 5, 6)
	    self.failUnlessEqual(imol_map, -1, "   bogus MTZ returns wrong molecule number")
	    now_n_molecules = graphics_n_molecules()
	    self.failUnlessEqual(now_n_molecules, pre_n_molecules, "   bogus MTZ creates extra map %s %s" %(pre_n_molecules, now_n_molecules))
	    
	    # correct mtz test
	    imol_map = make_and_draw_map(rnase_mtz, "FWT","PHWT","",0,0)
	    change_contour_level(0)
	    change_contour_level(0)
	    change_contour_level(0)
	    set_imol_refinement_map(imol_map)
	    self.failUnless(valid_map_molecule_qm(imol_map))

    def test07_0(self):
	    """Another Level Test"""
	    imol_map_2 = another_level()
	    self.failUnless(valid_map_molecule_qm(imol_map_2))
	    
    def test08_0(self):
	    """Set Atom Atribute Test"""
	    atom_ls = []
	    global imol_rnase
	    print "BL DEBUG:: set attr", imol_rnase,"A",11,""," CA ","","x",64.5
	    set_atom_attribute(imol_rnase,"A",11,""," CA ","","x",64.5) # an Angstrom or so
	    atom_ls = residue_info(imol_rnase, "A", 11, "")
	    print "BL DEBUG:: atom_ls", atom_ls
	    self.failIfEqual(atom_ls, [])
	    atom = atom_ls[0]
	    compound_name = atom[0]
	    atom_name = compound_name[0]
	    if (atom_name == " CA "):
		    x = atom_ls[0][2][0]
		    self.failUnlessAlmostEqual(x, 64.5)

    def test09_0(self):
	"""Add Terminal Residue Test"""
	import os, types
	#from types import *

	self.skipIf(not type(terminal_residue_test_pdb) is StringType, 
		"%s does not exist - skipping test" %terminal_residue_test_pdb)
	self.skipIf(not os.path.isfile(terminal_residue_test_pdb),
		"%s does not exist - skipping test" %terminal_residue_test_pdb)
	# OK, file exists
	imol = read_pdb(terminal_residue_test_pdb)
	self.skipIf(not valid_model_molecule_qm(imol), 
		"%s bad pdb read - skipping test" %terminal_residue_test_pdb)
	add_terminal_residue(imol, "A", 1, "ALA", 1)
	write_pdb_file(imol, "regression-test-terminal-residue.pdb")
	# where did that extra residue go?
	# 
	# Let's copy a fragment from imol and go to
	# the centre of that fragment, and check where
	# we are and compare it to where we expected
	# to be.
	# 
	new_mol = new_molecule_by_atom_selection(imol, "//A/0")
	print "BL DEBUG:: new mol ", new_mol
	self.failUnless(valid_model_molecule_qm(new_mol))
	move_molecule_here(new_mol)
	# not sure we want to move the molecule to the screen centre?!
	#set_go_to_atom_molecule(new_mol)
	#set_go_to_atom_chain_residue_atom_name("A", 0, " CA ")
	rc = rotation_centre()
	ls = [45.6, 15.8, 11.8]
	r = sum([rc[i] - ls[i] for i in range(len(rc))])
	print "BL DEBUG:: rc and r is ", rc, r
	self.failUnless(r>0.66, "Bad placement of terminal residue")

    def test10_0(self):
	    """Select by Sphere"""

	    global imol_rnase
	    imol_sphere = new_molecule_by_sphere_selection(imol_rnase, 
							   24.6114959716797, 24.8355808258057, 7.43978214263916,
							   3.6)

	    self.failUnless(valid_model_molecule_qm(imol_sphere), "Bad sphere molecule")
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
			

    def test11_0(self):
	"""Test Views"""
	view_number = add_view([32.0488, 21.6531, 13.7343],
			       [-0.12784, -0.491866, -0.702983, -0.497535],
			       20.3661,
			       "B11 View")
	go_to_view_number(view_number, 1)
	# test for something??

    def test12_0(self):
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

    def test13_0(self):
	    """Rotamer outliers"""

	    # pre-setup so that residues 1 and 2 are not rotamers "t" but 3
	    # to 8 are (1 and 2 are more than 40 degrees away). 9-16 are
	    # other residues and rotamers.
	    imol_rotamers = read_pdb(os.path.join(unittest_data_dir, "rotamer-test-fragment.pdb"))
				     
	    rotamer_anal = rotamer_graphs(imol_rotamers)
	    self.failUnless(type(rotamer_anal) is ListType)
	    self.failUnlessEqual(len(rotamer_anal), 14)
	    a_1 = rotamer_anal[0]
	    a_2 = rotamer_anal[1]
	    a_last = rotamer_anal[len(rotamer_anal)-1]

	    anal_str_a1 = a_1[-1]
	    anal_str_a2 = a_2[-1]
	    anal_str_a3 = a_last[-1]

	    self.failUnless((anal_str_a1 == "Rotamer not recognised" and
			     anal_str_a2 == "Rotamer not recognised" and
			     anal_str_a3 == "Missing Atoms"),
			    "  failure rotamer test: %s %s %s" %(a_1, a_2, a_last))
 
    # Don't reset the occupancies of the other parts of the residue
    # on regularization using alt confs
    #
    def test14_0(self):
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

	    self.failUnless(valid_model_molecule_qm(imol_fragment), "bad molecule for reading coords in Alt Conf Occ test")
	    occ_sum_pre = get_occ_sum(imol_fragment)

	    replace_state = refinement_immediate_replacement_state()
	    set_refinement_immediate_replacement(1)
	    regularize_zone(imol_fragment, "X", 15, 15, "A")
	    accept_regularizement()
	    if (replace_state == 0):
		    set_refinement_immediate_replacement(0)

	    occ_sum_post = get_occ_sum(imol_fragment)
	    self.failUnlessAlmostEqual(occ_sum_pre, occ_sum_post, 1, "   test for closeness: %s %s" %(occ_sum_pre, occ_sum_post))

    # This test we expect to fail until the CISPEP correction code is in
    # place (using mmdb-1.10+).
    # 
    def test15_0(self):
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

	    cis_pep_mol = read_pdb(os.path.join(unittest_data_dir, "tutorial-modern-cis-pep-12A_refmac0.pdb"))

	    view_number = add_view([63.455, 11.764, 1.268],
				   [-0.760536, -0.0910907, 0.118259, 0.631906],
				   15.7374,
				   "CIS-TRANS cispep View")
	    go_to_view_number(view_number, 1)

	    cis_trans_convert(cis_pep_mol, chain_id, resno, ins_code)
	    pepflip(cis_pep_mol, chain_id, resno, ins_code)
	    res_type = residue_name(cis_pep_mol, chain_id, resno, ins_code)
	    print "BL DEBUG:: res type", res_type
	    self.failUnless(type(res_type) is StringType)
	    mutate(cis_pep_mol, chain_id, resno, "", "GLY")
	    with_auto_accept(
		    [refine_zone, cis_pep_mol, chain_id, resno, (resno + 1), ""],
		    [accept_regularizement],
		    [mutate, cis_pep_mol, chain_id, resno, "", res_type],
		    [auto_fit_best_rotamer, resno, "", ins_code, chain_id, cis_pep_mol,
		     imol_refinement_map(), 1, 1],
		    [refine_zone, cis_pep_mol, chain_id, resno, (resno + 1), ""],
		    [accept_regularizement]
		    )
	    
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
	    print "BL DEBUG:: o list is", o
	    self.failUnlessEqual(len(o), 3)
	    print "CISPEPs: ", o
	    self.failUnlessEqual("3", parts)
	    
    def test16_0(self):
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
	    
    #def test17_0(self):
	#    """Setting multiple atom attributes"""
	    
