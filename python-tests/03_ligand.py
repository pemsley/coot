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

	imol = get_monomer("3GP")
	if (valid_model_molecule_qm(imol)):
            delete_residue_hydrogens(imol, "A", 1, "", "")
        self.failIf(not valid_model_molecule_qm(imol), "not valid molecule for 3GP")


    def test03_0(self):
        """Delete all-molecule Hydrogens"""
        imol = unittest_pdb("monomer-3GP.pdb")
        self.failUnless(valid_model_molecule_qm(imol))
        print "here 1"
        n = delete_hydrogens(imol)
        self.failUnless(n > 0)


    def test04_0(self):
        """Non-Autoloads"""

        def get_ccp4_version():
            s = shell_command_to_string("cad -i").split("\n")
            if s:
                for line in s:
                    if "CCP4 software suite: patch level" in line:
                        sl = line.split()
                        return sl[-1]
            return False

        def old_ccp4_restraints_qm():
            global have_ccp4_qm
            if (not have_ccp4_qm):
                return False
            else:
                v = get_ccp4_version()
                # will always be string
                return v < "6.2"

        r_1 = monomer_restraints("LIG")
        o = old_ccp4_restraints_qm()
        unittest_pdb("test-LIG.pdb")
        r_2 = monomer_restraints("LIG")
        remove_non_auto_load_residue_name("LIG")
        unittest_pdb("test-LIG.pdb")
        r_3 = monomer_restraints("LIG")
        delete_restraints("LIG")
        add_non_auto_load_residue_name("LIG")
        r_4 = monomer_restraints("LIG")
        # r_1, r_2, r_4 should be False, r_3 should be filled
        #
        # just to be clear here: either this is modern
        # restraints and we should get a list in r-3... or
        # this is old restraints. Either of those should
        # result in a PASS.
        self.failUnless(all([r_1 == False,
                             r_2 == False,
                             r_4 == False,
                             old_ccp4_restraints_qm() if \
                             old_ccp4_restraints_qm() else \
                             isinstance(r_3, dict)]))


    def test05_0(self):
        """Merge molecules of a ligand with a spec"""

        imol = unittest_pdb("tutorial-modern.pdb")
        imol_lig = get_monomer("3GP")

        spec = ["L", 1, ""]

        set_merge_molecules_ligand_spec(spec)
        merge_molecules([imol_lig], imol)

        # now check that L1 exists in imol

        # residue_spec_to_residue_name expects a 4-ele spec. Hmm
        spec[0:0] = [True]
        rn = residue_spec_to_residue_name(imol, spec)
        self.failUnless(isinstance(rn, str))
        self.failUnless(rn == "3GP")


    def test05_1(self):
        """Move and Refine Ligand test"""

        imol = unittest_pdb("tutorial-modern.pdb")
        self.failUnless(valid_model_molecule_qm(imol))
        new_rc = [55.3, 9.1, 20.6]
        # set the view
        view_number = add_view([54.5698, 8.7148, 20.5308],
                               [0.046229, -0.157139, -0.805581, 0.569395],
                               19.8858,
                               "ligand-view")

        go_to_view_number(view_number, 1)

        # update the map
        set_rotation_centre(*new_rc)
        move_molecule_here(imol)
        backup_mode = backup_state(imol)
        alt_conf = ""
        replacement_state = refinement_immediate_replacement_state()

        turn_off_backup(imol)
        set_refinement_immediate_replacement(1)
        refine_zone(imol, "A", 1, 1, alt_conf)
        accept_regularizement()
        rotate_y_scene(600, 0.1)
        if replacement_state == 0:
            set_refinement_immediate_replacement(0)
        if backup_mode == 1:
            turn_on_backup(imol)
        # ok. no fail.


    def test06_0(self):
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
        imol_map = auto_read_make_and_draw_maps(pdb43ca_mtz)
        imol_map_1 = imol_map[0]
        imol_map_2 = imol_map[1]

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


    def test07_0(self):
        """flip residue (around eigen vectors)"""

        # new version
        imol_orig = unittest_pdb("monomer-3GP.pdb")
        imol_copy = copy_molecule(imol_orig)

        self.failIf(not valid_model_molecule_qm(imol_orig),
                    "not valid molecule for monomer-3GP.pdb")

        # we need this, otherwise active-atom is (accidentally) the wrong
        # molecule
        set_go_to_atom_molecule(imol_copy)
        set_go_to_atom_chain_residue_atom_name("A", 1, " C8 ")

        active_atom = active_residue()
        self.failUnless(active_atom, "No active atom")
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        self.failIf(imol == imol_orig,
                    "oops - didn't pick the copy for active res")
        flip_ligand(imol, chain_id, res_no)
        atom_orig_1 = get_atom(imol_orig, "A", 1, "", " C8 ")
        atom_move_1 = get_atom(imol     , "A", 1, "", " C8 ")

        self.failUnless(isinstance(atom_orig_1, list), "atom_orig_1 not found")

        self.failUnless(isinstance(atom_move_1, list), "atom_move_1 not found")

        d = bond_length(atom_orig_1[2], atom_move_1[2])
        print "distance: ", d
        self.failUnless(d > 2.1, "fail to move test atom d1")
        flip_ligand(imol, chain_id, res_no)
        flip_ligand(imol, chain_id, res_no)
        flip_ligand(imol, chain_id, res_no)
        # having flipped it round the axes 4
        # times, we should be back where we
        # started.
        atom_orig_1 = get_atom(imol_orig, "A", 1, "", " C8 ")
        atom_move_1 = get_atom(imol     , "A", 1, "", " C8 ")
        d2 = bond_length(atom_orig_1[2], atom_move_1[2])
        print "distance d2: ", d2
        self.failUnless(d2 < 0.001, "fail to move atom back to start d2")


    def test08_0(self):
        """Test dipole"""

        if self.skip_test(True, "Skipping dipole test. disabled for now!"):
            return

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


    def test09_0(self):
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


    def test10_0(self):
        """Pyrogen Runs OK?"""

        if self.skip_test(not enhanced_ligand_coot_p(),
                          "No ligand enhanced version, skipping Pyrogen test"):
            return

        # bad things may well happen if we run the wrong version of pyrogen.
        # so force pyrogen to be the one that is installed alongside this version of coot
        # that we are running. We do that by looking and manipulating sys.argv[0]
        import os, sys

        coot_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
        prefix_dir = os.path.normpath(os.path.join(coot_dir, ".."))
        pyrogen_exe = "pyrogen"
        if is_windows():
            pyrogen_exe = "pyrogen.bat"
        pyrogen_bin = os.path.normpath(os.path.join(prefix_dir, "bin",
                                                    pyrogen_exe))

        smiles = "C1CNC1"
        tlc_text = "XXX"
        log_file_name = "pyrogen.log"

        # do we have pass now?
        if not enhanced_ligand_coot_p():
            # dont test pyrogen
            # not needed any more. Tested above.
            return
        else:
            global use_mogul

            if use_mogul:
                arg_list = ["--residue-type", tlc_text, smiles]
            else:
                if command_in_path_qm("mogul"):
                    arg_list = ["--residue-type", tlc_text, smiles]
                else:
                    arg_list = ["--no-mogul", "--residue-type", tlc_text, smiles]
                    # arg_list = ["--no-mogul", "-M", "--residue-type", tlc_text, smiles]
            popen_status = popen_command(pyrogen_bin, arg_list, [], log_file_name, True)
            # self.assertTrue(popen_status == 0)
            self.assertEqual(popen_status, 0,
                             "WARNING:: pyrogen exited with status %i\n" %popen_status)
            pdb_file_name = tlc_text + "-pyrogen.pdb"
            cif_file_name = tlc_text + "-pyrogen.cif"
            imol = handle_read_draw_molecule_with_recentre(pdb_file_name, 0)
            print "INFO:: pyrogen will try to read pdb file %s" %pdb_file_name
            # add test for chirality in the dictionary here
            self.assertTrue(valid_model_molecule_qm(imol))

    def test11_0(self):
        """pyrogen dictionary does not make double-quoted atom names"""

        # untested - no mogul

        # make sure that you are running the correct pyrogen

        if self.skip_test(not enhanced_ligand_coot_p(),
                          "No ligand enhanced version, skipping Pyrogen test"):
            return

        import os
        if os.path.isfile("UVP-pyrogen.cif"):
            os.remove("UVP-pyrogen.cif")

        popen_status = popen_command("pyrogen",
#                                     ["-nM", "-r", "UVP",
                                     ["-n", "-r", "UVP",
                                      "CO[C@@H]1[C@H](O)[C@H](O[C@H]1[n+]1ccc(O)nc1O)\\C=C\\P(O)(O)=O"],
                                     [], "pyrogen.log", False)
        self.assertEqual(popen_status, 0,
                         "Fail to correctly run pyrogen\n")
        read_cif_dictionary("UVP-pyrogen.cif")
        imol = get_monomer("UVP")
        self.failUnless(valid_model_molecule_qm(imol),
                        "Fail to load molecule from pyrogen dictionary\n")
        atom_info = residue_info(imol, "A", 1, "")
        # ok
        for atom in atom_info:
            atom_name = residue_atom2atom_name(atom)
            self.failIf("\"" in atom_name,
                        "Atom name quote fail %s" %atom_name)


    # FLEV will not make a PNG if it is not compiled with
    # C++-11 - and that is OK for 0.8.9.x.
    #
    # def test12_0(self):
    #     """FLEV makes a PNG"""

    #     import os

    #     if self.skip_test(not enhanced_ligand_coot_p(),
    #                       "No ligand enhanced version, skipping FLEV test"):
    #         return

    #     fn = "test-flev-greg-testcase.png"

    #     if os.path.exists(fn):
    #         os.remove(fn)

    #     imol = unittest_data_dir("tutorial-modern.pdb")
    #     imol_ligand = get_monomer("3GP")

    #     set_rotation_centre(54, 10, 20)
    #     move_molecule_to_screen_centre(imol_ligand)
    #     set_merge_molecules_ligand_spec(["L", 1, ""])
    #     merge_molecules([imol_ligand], imol)
    #     fle_view_with_rdkit_to_png(imol, "L", 1, "", 4.8, fn)
    #     self.failUnless(os.path.exists(fn))
