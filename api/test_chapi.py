
# test_molecules_container.py - Python-side tests for molecules_container_t
#
# These tests exercise the coot_headless_api nanobind module, mirroring
# selected tests from test-molecules-container.cc
#
# Usage:
#   MOORHEN_TEST_DATA_DIR=/path/to/test/data python -m unittest test_molecules_container
#

import coot_headless_api
import unittest
import json
import os

chapi = coot_headless_api.molecules_container_t(True)

def reference_data(filename):
    """Resolve test data file, honouring MOORHEN_TEST_DATA_DIR."""
    d = os.environ.get("MOORHEN_TEST_DATA_DIR", "")
    if d:
        return os.path.join(d, filename)
    return filename


@unittest.skipIf(coot_headless_api is None, "coot_headless_api not available")
class TestMoleculesContainer(unittest.TestCase):
    """Tests that exercise molecules_container_t via the nanobind API."""

    @classmethod
    def setUpClass(cls):
        chapi.geometry_init_standard()

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    def _read_tutorial_model(self, n=1):
        """Read moorhen tutorial structure n, return imol."""
        fn = reference_data(f"moorhen-tutorial-structure-number-{n}.pdb")
        imol = chapi.read_pdb(fn)
        self.assertTrue(chapi.is_valid_model_molecule(imol),
                        f"Failed to read {fn}")
        return imol

    def _read_tutorial_map(self, n=1):
        """Read moorhen tutorial map n, return imol_map."""
        fn = reference_data(f"moorhen-tutorial-map-number-{n}.mtz")
        imol_map = chapi.read_mtz(fn, "FWT", "PHWT", "W", False, False)
        self.assertTrue(chapi.is_valid_map_molecule(imol_map),
                        f"Failed to read {fn}")
        return imol_map

    # ------------------------------------------------------------------
    # test: read a coordinate file and verify basic properties
    # ------------------------------------------------------------------
    def test_read_coordinates(self):
        """Read a PDB file and check it has atoms and chains."""
        imol = self._read_tutorial_model(1)

        n_atoms = chapi.get_number_of_atoms(imol)
        self.assertGreater(n_atoms, 0, "Model should contain atoms")

        chains = chapi.get_chains_in_model(imol)
        self.assertGreater(len(chains), 0, "Model should contain chains")

        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: Ramachandran analysis
    #       (mirrors test_ramachandran_analysis in C++)
    # ------------------------------------------------------------------
    def test_ramachandran_analysis(self):
        """Ramachandran analysis should flag >100 favoured residues."""
        imol = self._read_tutorial_model(4)

        ra = chapi.ramachandran_analysis(imol)
        n_favoured = 0
        for chain_val in ra.cviv:
            for res_val in chain_val.rviv:
                if res_val.function_value > 0.5:
                    n_favoured += 1

        self.assertGreater(n_favoured, 100,
                           f"Expected >100 favoured residues, got {n_favoured}")
        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: mutate a residue
    #       (mirrors test_mutate in C++)
    # ------------------------------------------------------------------
    def test_mutate(self):
        """Mutate SER 270 to TYR and check the mutation succeeded."""
        imol = self._read_tutorial_model(4)

        cid = "//A/270/CA"
        mutate_status = chapi.mutate(imol, cid, "TYR")
        self.assertEqual(mutate_status, 1, "mutate() should return 1 on success")

        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: delete a residue using CID
    #       (mirrors test_weird_delete in C++)
    # ------------------------------------------------------------------
    def test_delete_using_cid(self):
        """Delete residue A/151 by CID and verify it is gone."""
        imol = self._read_tutorial_model(4)

        n_atoms_before = chapi.get_number_of_atoms(imol)
        chapi.delete_using_cid(imol, "//A/151", "RESIDUE")
        n_atoms_after = chapi.get_number_of_atoms(imol)

        self.assertLess(n_atoms_after, n_atoms_before,
                        "Atom count should decrease after deleting a residue")

        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: copy a fragment using a CID selection
    #       (mirrors test_copy_fragment_using_cid in C++)
    # ------------------------------------------------------------------
    def test_copy_fragment_using_cid(self):
        """Copy residues A/131-140 into a new molecule."""
        imol = self._read_tutorial_model(1)

        imol_new = chapi.copy_fragment_using_cid(imol, "//A/131-140")
        self.assertTrue(chapi.is_valid_model_molecule(imol_new),
                        "Copied fragment should be a valid molecule")

        n_atoms = chapi.get_number_of_atoms(imol_new)
        self.assertGreater(n_atoms, 0,
                           "Fragment molecule should contain atoms")

        chapi.close_molecule(imol_new)
        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: density correlation analysis
    #       (uses both model and map)
    # ------------------------------------------------------------------
    def test_density_correlation_analysis(self):
        """Per-residue density correlation should return results."""
        imol = self._read_tutorial_model(1)
        imol_map = self._read_tutorial_map(1)

        dca = chapi.density_correlation_analysis(imol, imol_map)
        n_residues = 0
        for chain_val in dca.cviv:
            n_residues += len(chain_val.rviv)

        self.assertGreater(n_residues, 100,
                           "Should have density correlation for >100 residues")

        chapi.close_molecule(imol_map)
        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: hydrogen atom manipulation
    # ------------------------------------------------------------------
    def test_add_and_delete_hydrogens(self):
        """Add hydrogens, verify count increases, then remove them."""
        imol = self._read_tutorial_model(1)

        n_h_before = chapi.get_number_of_hydrogen_atoms(imol)
        chapi.add_hydrogen_atoms(imol)
        n_h_after = chapi.get_number_of_hydrogen_atoms(imol)
        self.assertGreater(n_h_after, n_h_before,
                           "Hydrogen count should increase after add_hydrogen_atoms")

        chapi.delete_hydrogen_atoms(imol)
        n_h_final = chapi.get_number_of_hydrogen_atoms(imol)
        self.assertEqual(n_h_final, 0,
                         "Hydrogen count should be 0 after delete_hydrogen_atoms")

        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: molecule management (open/close/count)
    # ------------------------------------------------------------------
    def test_molecule_management(self):
        """Opening and closing molecules updates the count correctly."""
        n_start = chapi.get_number_of_molecules()

        imol1 = self._read_tutorial_model(1)
        imol2 = self._read_tutorial_model(4)
        self.assertNotEqual(imol1, imol2)

        chapi.close_molecule(imol1)
        chapi.close_molecule(imol2)

        self.assertFalse(chapi.is_valid_model_molecule(imol1))
        self.assertFalse(chapi.is_valid_model_molecule(imol2))

    # ------------------------------------------------------------------
    # test: rotamer analysis
    #       (mirrors test_rotamer_validation in C++)
    # ------------------------------------------------------------------
    def test_rotamer_analysis(self):
        """Rotamer analysis should find >150 residues with good rotamers."""
        imol = self._read_tutorial_model(1)

        ra = chapi.rotamer_analysis(imol)
        n_good = 0
        for chain_val in ra.cviv:
            for res_val in chain_val.rviv:
                if res_val.function_value > 0.5:
                    n_good += 1

        self.assertGreater(n_good, 150,
                           f"Expected >150 good rotamers, got {n_good}")
        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: density correlation values are sane
    #       (mirrors test_density_correlation_validation in C++)
    # ------------------------------------------------------------------
    def test_density_correlation_values_sane(self):
        """Density correlations should be in [-1, 1] with >250 good residues."""
        imol = self._read_tutorial_model(1)
        imol_map = self._read_tutorial_map(1)

        dca = chapi.density_correlation_analysis(imol, imol_map)
        n_good = 0
        bad_correls = False
        for chain_val in dca.cviv:
            for res_val in chain_val.rviv:
                if res_val.function_value > 0.5:
                    n_good += 1
                if res_val.function_value < -1.0 or res_val.function_value > 1.0:
                    bad_correls = True

        self.assertFalse(bad_correls, "No correlation values should be outside [-1, 1]")
        self.assertGreater(n_good, 250,
                           f"Expected >250 residues with good density correlation, got {n_good}")

        chapi.close_molecule(imol_map)
        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: auto-fit rotamer
    #       (mirrors test_auto_fit_rotamer_1 in C++)
    # ------------------------------------------------------------------
    def test_auto_fit_rotamer(self):
        """Auto-fit rotamer for A/61 should return success."""
        imol = self._read_tutorial_model(4)
        imol_map = self._read_tutorial_map(4)

        status = chapi.auto_fit_rotamer(imol, "A", 61, "", "", imol_map)
        self.assertEqual(status, 1, "auto_fit_rotamer should return 1 on success")

        chapi.close_molecule(imol_map)
        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: difference map peaks
    #       (mirrors test_difference_map_peaks in C++)
    # ------------------------------------------------------------------
    def test_difference_map_peaks(self):
        """Should find >8 difference map peaks at 4.5 sigma."""
        imol = self._read_tutorial_model(1)
        imol_diff_map = chapi.read_mtz(
            reference_data("moorhen-tutorial-map-number-1.mtz"),
            "DELFWT", "PHDELWT", "W", False, True)
        self.assertTrue(chapi.is_valid_map_molecule(imol_diff_map))

        n_rmsd = 4.5
        sites = chapi.difference_map_peaks(imol_diff_map, imol, n_rmsd)
        self.assertGreater(len(sites), 8,
                           f"Expected >8 difference map peaks, got {len(sites)}")

        chapi.close_molecule(imol_diff_map)
        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: get monomer
    # ------------------------------------------------------------------
    def test_get_monomer(self):
        """Get a monomer from the dictionary and verify it is valid."""
        imol = chapi.get_monomer("TYR")
        self.assertTrue(chapi.is_valid_model_molecule(imol),
                        "get_monomer('TYR') should return a valid molecule")

        n_atoms = chapi.get_number_of_atoms(imol)
        self.assertGreater(n_atoms, 0, "TYR monomer should have atoms")

        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: add waters
    #       (mirrors test_add_water in C++)
    # ------------------------------------------------------------------
    def test_add_waters(self):
        """Adding waters to a model should find >10 waters."""
        imol = self._read_tutorial_model(1)
        imol_map = self._read_tutorial_map(1)

        n_waters = chapi.add_waters(imol, imol_map)
        self.assertGreater(n_waters, 10,
                           f"Expected >10 waters to be added, got {n_waters}")

        chapi.close_molecule(imol_map)
        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: get single-letter codes for chain
    # ------------------------------------------------------------------
    def test_get_single_letter_codes(self):
        """Sequence extraction should return a non-empty string."""
        imol = self._read_tutorial_model(1)

        chains = chapi.get_chains_in_model(imol)
        self.assertGreater(len(chains), 0)

        seq = chapi.get_single_letter_codes_for_chain(imol, chains[0])
        self.assertGreater(len(seq), 50,
                           "Should get a substantial sequence string")

        chapi.close_molecule(imol)


    # ------------------------------------------------------------------
    # test: flip peptide and undo
    #       (mirrors test_undo_and_redo in C++)
    # ------------------------------------------------------------------
    def test_flip_peptide_and_undo(self):
        """Flip a peptide, undo it, atom count should be preserved."""
        imol = self._read_tutorial_model(1)

        n_atoms_before = chapi.get_number_of_atoms(imol)
        chapi.flip_peptide_using_cid(imol, "//A/14/CA", "")
        n_atoms_after_flip = chapi.get_number_of_atoms(imol)
        self.assertEqual(n_atoms_before, n_atoms_after_flip,
                         "Peptide flip should not change atom count")

        chapi.undo(imol)
        n_atoms_after_undo = chapi.get_number_of_atoms(imol)
        self.assertEqual(n_atoms_before, n_atoms_after_undo,
                         "Undo should preserve atom count")

        chapi.close_molecule(imol)

    # ------------------------------------------------------------------
    # test: reading a missing file should fail gracefully
    #       (mirrors test_read_a_missing_map in C++)
    # ------------------------------------------------------------------
    def test_read_missing_file(self):
        """Reading a non-existent file should return an invalid molecule."""
        imol = chapi.read_pdb("a-file-that-does-not-exist.pdb")
        self.assertFalse(chapi.is_valid_model_molecule(imol))

        imol_map = chapi.read_ccp4_map("a-map-that-does-not-exist.map", False)
        self.assertFalse(chapi.is_valid_map_molecule(imol_map))

    # ------------------------------------------------------------------
    # test: SSM superpose
    # ------------------------------------------------------------------
    def test_SSM_superpose(self):
        """Superposing a structure onto itself should give near-zero RMSD."""
        imol1 = self._read_tutorial_model(1)
        imol2 = self._read_tutorial_model(1)
        result = chapi.SSM_superpose(imol1, "A", imol2, "A")
        results_dict = json.loads(result.superpose_info)
        rmsd = results_dict["rmsd"]
        print("      rmsd", rmsd)

        # After superposition the RMSD should be very small (same structure)
        self.assertAlmostEqual(rmsd, 0.0, places=2, msg=f"Self-superposition RMSD should be ~0, got {rmsd}")

        chapi.close_molecule(imol2)
        chapi.close_molecule(imol1)


if __name__ == "__main__":
    unittest.main()
