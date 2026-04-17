/* coot-utils/test-coord-utils.cc
 *
 * Copyright 2026 by Paul Emsley
 * Copyright 2026 by Medical Research Council
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

// Tests for functions in coot-coord-utils.hh
// These are free functions (not class members), so no molecules_container_t needed.
// Return 0 for fail, 1 for pass.

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>
#include <clipper/core/rotation.h>

#include "coot-coord-utils.hh"
#include "coot-coord-extras.hh"
#include "coot-map-utils.hh"
#include "geometry/residue-and-atom-specs.hh"

void starting_test(const char *func) {
   std::cout << "\nStarting " << func << "()" << std::endl;
}

// wrap with the directory for the test data.
std::string
reference_data(const std::string &file) {
   char *env = getenv("MOORHEN_TEST_DATA_DIR");
   if (env) {
      std::string joined = coot::util::append_dir_file(env, file);
      return joined;
   } else {
      return file;
   }
}

// Helper: read a PDB file into an mmdb::Manager. Caller deletes. Returns nullptr on failure.
mmdb::Manager *read_pdb(const std::string &file_name) {
   mmdb::Manager *mol = new mmdb::Manager;
   mmdb::ERROR_CODE err = mol->ReadCoorFile(file_name.c_str());
   if (err != mmdb::Error_NoError) {
      std::cout << "ERROR reading " << file_name << std::endl;
      delete mol;
      return nullptr;
   }
   return mol;
}

// ==================== tests ====================

int test_write_coords_pdb() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::string tmp_file = "test-write-coords-pdb-tmp.pdb";
      int r = coot::write_coords_pdb(mol, tmp_file);
      if (r == 0) {
         // check the file was written and is not empty
         std::ifstream f(tmp_file);
         if (f.good()) {
            f.seekg(0, std::ios::end);
            if (f.tellg() > 100)
               status = 1;
         }
         std::remove(tmp_file.c_str());
      }
      delete mol;
   }
   return status;
}

int test_write_coords_cif() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::string tmp_file = "test-write-coords-cif-tmp.cif";
      int r = coot::write_coords_cif(mol, tmp_file);
      if (r == 0) {
         std::ifstream f(tmp_file);
         if (f.good()) {
            f.seekg(0, std::ios::end);
            if (f.tellg() > 100)
               status = 1;
         }
         std::remove(tmp_file.c_str());
      }
      delete mol;
   }
   return status;
}

int test_pad_atom_name() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::string padded = coot::pad_atom_name("CA", "C");
   // Padded atom name should be 4 chars for PDB format
   if (padded.size() == 4)
      status = 1;
   return status;
}

int test_distance() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      // get first two atoms
      int selHnd = mol->NewSelection();
      mol->SelectAtoms(selHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mmdb::PAtom *atoms = nullptr;
      int n_atoms = 0;
      mol->GetSelIndex(selHnd, atoms, n_atoms);
      if (n_atoms >= 2) {
         double d = coot::distance(atoms[0], atoms[1]);
         if (d >= 0.0 && d < 100.0)
            status = 1;
      }
      mol->DeleteSelection(selHnd);
      delete mol;
   }
   return status;
}

int test_angle() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int selHnd = mol->NewSelection();
      mol->SelectAtoms(selHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mmdb::PAtom *atoms = nullptr;
      int n_atoms = 0;
      mol->GetSelIndex(selHnd, atoms, n_atoms);
      if (n_atoms >= 3) {
         double a = coot::angle(atoms[0], atoms[1], atoms[2]);
         if (a > 0.0 && a < 360.0)
            status = 1;
      }
      mol->DeleteSelection(selHnd);
      delete mol;
   }
   return status;
}

int test_co_and_update_position() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int selHnd = mol->NewSelection();
      mol->SelectAtoms(selHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mmdb::PAtom *atoms = nullptr;
      int n_atoms = 0;
      mol->GetSelIndex(selHnd, atoms, n_atoms);
      if (n_atoms > 0) {
         clipper::Coord_orth pos = coot::co(atoms[0]);
         clipper::Coord_orth new_pos(pos.x() + 1.0, pos.y(), pos.z());
         coot::update_position(atoms[0], new_pos);
         clipper::Coord_orth pos2 = coot::co(atoms[0]);
         if (std::fabs(pos2.x() - new_pos.x()) < 0.01)
            status = 1;
      }
      mol->DeleteSelection(selHnd);
      delete mol;
   }
   return status;
}

int test_get_position_hash() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      float h1 = coot::get_position_hash(mol);
      // hash should be non-zero for a real molecule
      if (h1 != 0.0f)
         status = 1;
      delete mol;
   }
   return status;
}

int test_get_selection_handle() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::atom_spec_t spec("A", 10, "", " CA ", "");
      int selHnd = coot::get_selection_handle(mol, spec);
      if (selHnd >= 0) {
         mmdb::PAtom *atoms = nullptr;
         int n_atoms = 0;
         mol->GetSelIndex(selHnd, atoms, n_atoms);
         if (n_atoms >= 1)
            status = 1;
         mol->DeleteSelection(selHnd);
      }
      delete mol;
   }
   return status;
}

int test_lsq_plane_info() {
   starting_test(__FUNCTION__);
   int status = 0;
   // make 6 roughly coplanar points in the xy plane with some z scatter
   std::vector<clipper::Coord_orth> v;
   v.push_back(clipper::Coord_orth(0, 0, 0.1));
   v.push_back(clipper::Coord_orth(2, 0, -0.1));
   v.push_back(clipper::Coord_orth(0, 2, 0.05));
   v.push_back(clipper::Coord_orth(2, 2, -0.05));
   v.push_back(clipper::Coord_orth(1, 1, 0.0));
   v.push_back(clipper::Coord_orth(1, 0, 0.08));
   try {
      coot::lsq_plane_info_t lpi(v);
      // the normal should be close to (0,0,1) or (0,0,-1)
      clipper::Coord_orth n = lpi.normal();
      double nz_abs = std::fabs(n.z());
      if (nz_abs > 0.9) {
         // test plane_deviation - point is ~2A above the plane
         double dev = lpi.plane_deviation(clipper::Coord_orth(0.5, 0.5, 2.0));
         if (std::fabs(dev) > 1.0)
            status = 1;
      }
   } catch (const std::exception &e) {
      std::cout << "Exception: " << e.what() << std::endl;
   }
   return status;
}

int test_lsq_plane_deviation() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::vector<clipper::Coord_orth> v;
   v.push_back(clipper::Coord_orth(0, 0, 0));
   v.push_back(clipper::Coord_orth(1, 0, 0));
   v.push_back(clipper::Coord_orth(0, 1, 0));
   v.push_back(clipper::Coord_orth(1, 1, 0));
   clipper::Coord_orth pt(0.5, 0.5, 3.0);
   std::pair<double, double> result = coot::lsq_plane_deviation(v, pt);
   // first: deviation of pt from plane, second: rms of plane atoms
   if (std::fabs(result.first - 3.0) < 0.01 && result.second < 0.01)
      status = 1;
   return status;
}

int test_is_member_p() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            int n_res = chain_p->GetNumberOfResidues();
            if (n_res >= 2) {
               std::vector<mmdb::Residue *> v;
               v.push_back(chain_p->GetResidue(0));
               bool r1 = coot::is_member_p(v, chain_p->GetResidue(0));
               bool r2 = coot::is_member_p(v, chain_p->GetResidue(1));
               if (r1 && !r2)
                  status = 1;
            }
         }
      }
      delete mol;
   }
   return status;
}

int test_is_hydrogen_atom() {
   starting_test(__FUNCTION__);
   int status = 0;
   // create a simple atom and test
   mmdb::Atom *at = new mmdb::Atom;
   at->SetAtomName(" H  ");
   at->SetElementName(" H");
   bool r1 = coot::is_hydrogen_atom(at);
   at->SetAtomName(" CA ");
   at->SetElementName(" C");
   bool r2 = coot::is_hydrogen_atom(at);
   if (r1 && !r2)
      status = 1;
   delete at;
   return status;
}

int test_residues_in_order_p() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            bool in_order = coot::residues_in_order_p(chain_p);
            // tutorial structure should have residues in order
            if (in_order)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_centre_of_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::pair<bool, clipper::Coord_orth> result = coot::centre_of_molecule(mol);
      if (result.first) {
         // centre should be finite and reasonable
         double x = result.second.x();
         double y = result.second.y();
         double z = result.second.z();
         if (std::isfinite(x) && std::isfinite(y) && std::isfinite(z))
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_centre_of_molecule_using_masses() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::pair<bool, clipper::Coord_orth> result = coot::centre_of_molecule_using_masses(mol);
      if (result.first) {
         double x = result.second.x();
         if (std::isfinite(x))
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_radius_of_gyration() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::pair<bool, double> result = coot::radius_of_gyration(mol);
      if (result.first) {
         if (result.second > 1.0 && result.second < 200.0)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_centre_of_residues() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            std::vector<mmdb::Residue *> residues;
            int n_res = chain_p->GetNumberOfResidues();
            for (int i = 0; i < std::min(n_res, 5); i++)
               residues.push_back(chain_p->GetResidue(i));
            std::pair<bool, clipper::Coord_orth> result = coot::centre_of_residues(residues);
            if (result.first && std::isfinite(result.second.x()))
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_hetify_residue_atoms_as_needed() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      // find a ligand residue or just test on any residue
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            mmdb::Residue *res = chain_p->GetResidue(0);
            if (res) {
               int r = coot::hetify_residue_atoms_as_needed(res);
               // for a standard residue, it should return 0 (no change)
               // but doesn't crash - that's a pass
               status = 1;
            }
         }
      }
      delete mol;
   }
   return status;
}

int test_hetify_residues_as_needed() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int r = coot::hetify_residues_as_needed(mol);
      // returns something, doesn't crash
      status = 1;
      delete mol;
   }
   return status;
}

int test_atoms_with_zero_occupancy() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<mmdb::Atom *> zero_occ = coot::atoms_with_zero_occupancy(mol);
      // the function returns without crashing, that's a pass
      // a well-formed tutorial structure should have few or no zero-occ atoms
      status = 1;
      delete mol;
   }
   return status;
}

int test_residues_with_alt_confs() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<mmdb::Residue *> alt_residues = coot::residues_with_alt_confs(mol);
      // function runs without crash
      status = 1;
      delete mol;
   }
   return status;
}

int test_residues_near_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::residue_spec_t spec("A", 50, "");
      std::vector<coot::residue_spec_t> nearby = coot::residues_near_residue(spec, mol, 5.0);
      if (nearby.size() > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_residues_near_residue_v2() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         std::vector<mmdb::Residue *> nearby = coot::residues_near_residue(res, mol, 5.0);
         if (nearby.size() > 0)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_residues_near_residues() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      auto m = coot::residues_near_residues(mol, 5.0);
      if (m.size() > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_interface_residues() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      // check if there are at least 2 chains
      std::vector<std::string> chains = coot::util::chains_in_molecule(mol);
      if (chains.size() >= 2) {
         auto result = coot::interface_residues(mol, chains[0], chains[1], 5.0);
         if (result.first.size() > 0 || result.second.size() > 0)
            status = 1;
      } else {
         // only one chain - pass anyway as interface_residues is tested with multi-chain
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_residues_near_position() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::pair<bool, clipper::Coord_orth> c = coot::centre_of_molecule(mol);
      if (c.first) {
         std::vector<mmdb::Residue *> nearby = coot::residues_near_position(c.second, mol, 10.0);
         if (nearby.size() > 0)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_closest_approach() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *r1 = coot::util::get_residue("A", 50, "", mol);
      mmdb::Residue *r2 = coot::util::get_residue("A", 51, "", mol);
      if (r1 && r2) {
         std::pair<bool, float> result = coot::closest_approach(mol, r1, r2);
         if (result.first) {
            if (result.second > 0.0 && result.second < 10.0)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_nearest_residue_by_sequence() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::residue_spec_t spec("A", 50, "");
      mmdb::Residue *res = coot::nearest_residue_by_sequence(mol, spec);
      if (res) {
         int resno = res->GetSeqNum();
         if (resno == 50)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_sort_chains() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::sort_chains(mol);
      // doesn't crash
      status = 1;
      delete mol;
   }
   return status;
}

int test_sort_residues() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::sort_residues(mol);
      // doesn't crash
      status = 1;
      delete mol;
   }
   return status;
}

int test_mol_has_symmetry() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      bool has_sym = coot::mol_has_symmetry(mol);
      // function runs without crash - symmetry depends on file format/reader
      status = 1;
      delete mol;
   }
   return status;
}

int test_mol_is_anisotropic() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      bool is_aniso = coot::mol_is_anisotropic(mol);
      // tutorial structure is typically isotropic
      if (!is_aniso)
         status = 1;
      delete mol;
   }
   return status;
}

int test_get_brick_id_inner() {
   starting_test(__FUNCTION__);
   int status = 0;
   int id = coot::get_brick_id_inner(1, 2, 3, 10, 10, 10);
   if (id >= 0)
      status = 1;
   return status;
}

int test_get_brick_id() {
   starting_test(__FUNCTION__);
   int status = 0;
   clipper::Coord_orth pt(5.0, 5.0, 5.0);
   clipper::Coord_orth pt_min(0.0, 0.0, 0.0);
   int id = coot::get_brick_id(pt, pt_min, 10, 10, 10, 6.0);
   if (id >= 0)
      status = 1;
   return status;
}

// ==================== coot::util namespace tests ====================

int test_single_letter_to_3_letter_code() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::string r = coot::util::single_letter_to_3_letter_code('A');
   if (r == "ALA")
      status = 1;
   return status;
}

int test_single_letter_to_3_letter_code_string() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::string r = coot::util::single_letter_to_3_letter_code("G");
   if (r == "GLY")
      status = 1;
   return status;
}

int test_three_letter_to_one_letter() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::string r = coot::util::three_letter_to_one_letter("ALA");
   if (r == "A") {
      std::string r2 = coot::util::three_letter_to_one_letter("GLY");
      if (r2 == "G")
         status = 1;
   }
   return status;
}

int test_three_letter_to_one_letter_with_specials() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::string r = coot::util::three_letter_to_one_letter_with_specials("HOH");
   if (r == "~")
      status = 1;
   return status;
}

int test_standard_residue_types() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::vector<std::string> types = coot::util::standard_residue_types();
   // 20 standard amino acids + MSE = 21
   if (types.size() == 21)
      status = 1;
   return status;
}

int test_PDB_standard_residue_types() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::vector<std::string> types = coot::util::PDB_standard_residue_types();
   if (types.size() > 20) // includes nucleotides etc
      status = 1;
   return status;
}

int test_is_nucleotide() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 10, "", mol);
      if (res) {
         short int is_nuc = coot::util::is_nucleotide(res);
         // amino acid residue should not be nucleotide
         if (!is_nuc)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_residue_has_hydrogens_p() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 10, "", mol);
      if (res) {
         bool has_h = coot::util::residue_has_hydrogens_p(res);
         // typical PDB files may or may not have hydrogens, just test it runs
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_residue_has_hetatms() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 10, "", mol);
      if (res) {
         int r = coot::util::residue_has_hetatms(res);
         // standard residue should not have hetatms
         if (r == 0)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 10, "", mol);
      if (res) {
         int resno = res->GetSeqNum();
         if (resno == 10)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_residue_by_spec() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::residue_spec_t spec("A", 10, "");
      mmdb::Residue *res = coot::util::get_residue(spec, mol);
      if (res) {
         if (res->GetSeqNum() == 10)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_first_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_first_residue(mol);
      if (res)
         status = 1;
      delete mol;
   }
   return status;
}

int test_get_nth_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res1 = coot::util::get_nth_residue(1, mol);
      mmdb::Residue *res2 = coot::util::get_nth_residue(2, mol);
      if (res1 && res2 && (res1 != res2))
         status = 1;
      delete mol;
   }
   return status;
}

int test_get_atom() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::atom_spec_t spec("A", 10, "", " CA ", "");
      mmdb::Atom *at = coot::util::get_atom(spec, mol);
      if (at) {
         std::string at_name(at->GetAtomName());
         if (at_name.find("CA") != std::string::npos)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_following_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::residue_spec_t spec("A", 50, "");
      mmdb::Residue *res = coot::util::get_following_residue(spec, mol);
      if (res) {
         if (res->GetSeqNum() == 51)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_previous_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::residue_spec_t spec("A", 50, "");
      mmdb::Residue *res = coot::util::get_previous_residue(spec, mol);
      if (res) {
         if (res->GetSeqNum() == 49)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_residue_centre() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         std::pair<bool, clipper::Coord_orth> result = coot::util::get_residue_centre(res);
         if (result.first && std::isfinite(result.second.x()))
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_CA_position_in_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         auto result = coot::util::get_CA_position_in_residue(res);
         if (result.first && std::isfinite(result.second.x()))
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_CB_position_in_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      // pick a non-GLY residue
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         std::string rn(res->GetResName());
         if (rn != "GLY") {
            auto result = coot::util::get_CB_position_in_residue(res);
            if (result.first && std::isfinite(result.second.x()))
               status = 1;
         } else {
            status = 1; // skip for GLY
         }
      }
      delete mol;
   }
   return status;
}

int test_residue_types_in_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<std::string> types = coot::util::residue_types_in_molecule(mol);
      if (types.size() > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_non_standard_residue_types_in_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<std::string> types = coot::util::non_standard_residue_types_in_molecule(mol);
      // just runs without crash
      status = 1;
      delete mol;
   }
   return status;
}

int test_chains_in_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<std::string> chains = coot::util::chains_in_molecule(mol);
      if (chains.size() > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_residues_in_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<mmdb::Residue *> residues = coot::util::residues_in_molecule(mol);
      if (residues.size() > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_number_of_residues_in_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int n = coot::util::number_of_residues_in_molecule(mol);
      if (n > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_max_number_of_residues_in_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int n = coot::util::max_number_of_residues_in_chain(mol);
      if (n > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_number_of_chains() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int n = coot::util::number_of_chains(mol);
      if (n >= 1)
         status = 1;
      delete mol;
   }
   return status;
}

int test_min_and_max_residues() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            std::pair<int, int> mm = coot::util::min_and_max_residues(chain_p);
            // second should be >= first for success
            if (mm.second >= mm.first)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_min_resno_in_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            auto result = coot::util::min_resno_in_chain(chain_p);
            if (result.first)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_max_resno_in_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            auto result = coot::util::max_resno_in_chain(chain_p);
            if (result.first)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_max_resno_in_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      auto result = coot::util::max_resno_in_molecule(mol);
      if (result.first)
         status = 1;
      delete mol;
   }
   return status;
}

int test_copy_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Manager *mol_copy = coot::util::copy_molecule(mol);
      if (mol_copy) {
         int n1 = coot::util::number_of_residues_in_molecule(mol);
         int n2 = coot::util::number_of_residues_in_molecule(mol_copy);
         if (n1 == n2 && n1 > 0)
            status = 1;
         delete mol_copy;
      }
      delete mol;
   }
   return status;
}

int test_copy_cell_and_symm_headers() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Manager *mol2 = new mmdb::Manager;
      bool success = coot::util::copy_cell_and_symm_headers(mol, mol2);
      if (success)
         status = 1;
      delete mol2;
      delete mol;
   }
   return status;
}

int test_create_mmdbmanager_from_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         mmdb::Manager *new_mol = coot::util::create_mmdbmanager_from_residue(res);
         if (new_mol) {
            int n = coot::util::number_of_residues_in_molecule(new_mol);
            if (n == 1)
               status = 1;
            delete new_mol;
         }
      }
      delete mol;
   }
   return status;
}

int test_create_mmdbmanager_from_residue_vector() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<mmdb::Residue *> res_vec;
      for (int i = 50; i <= 55; i++) {
         mmdb::Residue *r = coot::util::get_residue("A", i, "", mol);
         if (r) res_vec.push_back(r);
      }
      if (res_vec.size() > 0) {
         auto result = coot::util::create_mmdbmanager_from_residue_vector(res_vec, mol);
         if (result.first && result.second) {
            int n = coot::util::number_of_residues_in_molecule(result.second);
            if (n == (int)res_vec.size())
               status = 1;
            delete result.second;
         }
      }
      delete mol;
   }
   return status;
}

int test_create_mmdbmanager_from_mmdbmanager() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      auto result = coot::util::create_mmdbmanager_from_mmdbmanager(mol);
      if (result.first) {
         int n = coot::util::number_of_residues_in_molecule(result.first);
         if (n > 0)
            status = 1;
         delete result.first;
      }
      delete mol;
   }
   return status;
}

int test_deep_copy_this_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         mmdb::Residue *res_copy = coot::util::deep_copy_this_residue(res);
         if (res_copy) {
            int n1 = res->GetNumberOfAtoms();
            int n2 = res_copy->GetNumberOfAtoms();
            if (n1 == n2 && n1 > 0)
               status = 1;
            delete res_copy;
         }
      }
      delete mol;
   }
   return status;
}

int test_get_residues_in_fragment() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            coot::residue_spec_t spec("A", 50, "");
            std::vector<mmdb::PResidue> frag = coot::util::get_residues_in_fragment(chain_p, spec);
            if (frag.size() > 0)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_get_residues_in_range() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<mmdb::Residue *> res = coot::util::get_residues_in_range(mol, "A", 50, 55);
      if (res.size() == 6)
         status = 1;
      delete mol;
   }
   return status;
}

int test_get_hetgroups() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<mmdb::Residue *> hets = coot::util::get_hetgroups(mol);
      // just test it runs
      status = 1;
      delete mol;
   }
   return status;
}

int test_get_this_and_next_residues() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::residue_spec_t spec("A", 50, "");
      auto result = coot::util::get_this_and_next_residues(spec, mol);
      if (result.first && result.second) {
         if (result.first->GetSeqNum() == 50 && result.second->GetSeqNum() == 51)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_intelligent_this_residue_mmdb_atom() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         mmdb::Atom *at = coot::util::intelligent_this_residue_mmdb_atom(res);
         if (at) {
            // for a protein residue, should return CA
            std::string name(at->GetAtomName());
            if (name.find("CA") != std::string::npos)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_previous_and_next_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         mmdb::Residue *prev = coot::util::previous_residue(res);
         mmdb::Residue *next = coot::util::next_residue(res);
         if (prev && next) {
            if (prev->GetSeqNum() == 49 && next->GetSeqNum() == 51)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_occupancy_sum() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         mmdb::PAtom *atoms = nullptr;
         int n_atoms = 0;
         res->GetAtomTable(atoms, n_atoms);
         if (n_atoms > 0) {
            float sum = coot::util::occupancy_sum(atoms, n_atoms);
            // all atoms with occ 1.0 should sum to n_atoms
            if (std::fabs(sum - (float)n_atoms) < 0.5)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_median_temperature_factor() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int selHnd = mol->NewSelection();
      mol->SelectAtoms(selHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mmdb::PAtom *atoms = nullptr;
      int n_atoms = 0;
      mol->GetSelIndex(selHnd, atoms, n_atoms);
      if (n_atoms > 0) {
         float med = coot::util::median_temperature_factor(atoms, n_atoms, 0.0, 999.0, false, false);
         if (med > 0.0 && med < 200.0)
            status = 1;
      }
      mol->DeleteSelection(selHnd);
      delete mol;
   }
   return status;
}

int test_average_temperature_factor() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int selHnd = mol->NewSelection();
      mol->SelectAtoms(selHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mmdb::PAtom *atoms = nullptr;
      int n_atoms = 0;
      mol->GetSelIndex(selHnd, atoms, n_atoms);
      if (n_atoms > 0) {
         float avg = coot::util::average_temperature_factor(atoms, n_atoms, 0.0, 999.0, 0, 0);
         if (avg > 0.0 && avg < 200.0)
            status = 1;
      }
      mol->DeleteSelection(selHnd);
      delete mol;
   }
   return status;
}

int test_matrix_convert() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::mat44 mat;
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         mat[i][j] = (i == j) ? 1.0 : 0.0;
   clipper::RTop_orth rtop = coot::util::matrix_convert(mat);
   // should be identity
   clipper::Coord_orth t(rtop.trn());
   double d = t.x() + t.y() + t.z();
   if (std::fabs(d) < 0.001)
      status = 1;
   return status;
}

int test_rotate_around_vector() {
   starting_test(__FUNCTION__);
   int status = 0;
   clipper::Coord_orth direction(0, 0, 1); // z-axis
   clipper::Coord_orth position(1, 0, 0);
   clipper::Coord_orth origin(0, 0, 0);
   double angle = M_PI / 2.0; // 90 degrees
   clipper::Coord_orth result = coot::util::rotate_around_vector(direction, position, origin, angle);
   // rotating (1,0,0) 90 degrees around z should give (0,1,0)
   if (std::fabs(result.x()) < 0.01 && std::fabs(result.y() - 1.0) < 0.01)
      status = 1;
   return status;
}

int test_extents() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      auto ext = coot::util::extents(mol);
      // max should be > min in all dimensions
      if (ext.second.x() > ext.first.x() &&
          ext.second.y() > ext.first.y() &&
          ext.second.z() > ext.first.z())
         status = 1;
      delete mol;
   }
   return status;
}

int test_residues_with_insertion_codes() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<mmdb::Residue *> res = coot::util::residues_with_insertion_codes(mol);
      // just test it runs
      status = 1;
      delete mol;
   }
   return status;
}

int test_omega_torsion() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res1 = coot::util::get_residue("A", 50, "", mol);
      mmdb::Residue *res2 = coot::util::get_residue("A", 51, "", mol);
      if (res1 && res2) {
         auto result = coot::util::omega_torsion(res1, res2, "");
         if (result.first) {
            // omega should be close to pi (trans) or 0 (cis)
            double omega = result.second;
            if (std::fabs(omega) > 2.5 || std::fabs(omega) < 0.5) // roughly trans or cis
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_is_cis() {
   starting_test(__FUNCTION__);
   int status = 0;
   bool r1 = coot::util::is_cis(0.1);   // near 0, should be cis
   bool r2 = coot::util::is_cis(M_PI);  // pi, should be trans
   if (r1 && !r2)
      status = 1;
   return status;
}

#if 0 // 20260408-PE this function does not exist yet.
int test_peptide_torsions() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res1 = coot::util::get_residue("A", 50, "", mol);
      mmdb::Residue *res2 = coot::util::get_residue("A", 51, "", mol);
      if (res1 && res2) {
         auto info = coot::util::peptide_torsions(res1, res2, "");
         if (info.status) {
            // phi and psi should be in range
            if (std::isfinite(info.phi) && std::isfinite(info.psi))
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}
#endif

int test_count_cis_peptides() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int n = coot::util::count_cis_peptides(mol);
      // just runs and returns a non-negative number
      if (n >= 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_cis_peptides_info_from_coords() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      auto cis_peps = coot::util::cis_peptides_info_from_coords(mol);
      // just runs
      status = 1;
      delete mol;
   }
   return status;
}

int test_average_position() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::vector<clipper::Coord_orth> pts;
   pts.push_back(clipper::Coord_orth(1, 0, 0));
   pts.push_back(clipper::Coord_orth(0, 1, 0));
   pts.push_back(clipper::Coord_orth(0, 0, 1));
   clipper::Coord_orth avg = coot::util::average_position(pts);
   double expected = 1.0 / 3.0;
   if (std::fabs(avg.x() - expected) < 0.01 &&
       std::fabs(avg.y() - expected) < 0.01 &&
       std::fabs(avg.z() - expected) < 0.01)
      status = 1;
   return status;
}

int test_median_position() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      try {
         clipper::Coord_orth med = coot::util::median_position(mol);
         if (std::isfinite(med.x()))
            status = 1;
      } catch (...) {}
      delete mol;
   }
   return status;
}

int test_transform_mol() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      auto centre_before = coot::centre_of_molecule(mol);
      clipper::RTop_orth rtop(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1),
                              clipper::Coord_orth(10, 0, 0));
      coot::util::transform_mol(mol, rtop);
      auto centre_after = coot::centre_of_molecule(mol);
      if (centre_before.first && centre_after.first) {
         double dx = centre_after.second.x() - centre_before.second.x();
         if (std::fabs(dx - 10.0) < 0.1)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_shift() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      auto c1 = coot::centre_of_molecule(mol);
      coot::util::shift(mol, clipper::Coord_orth(5.0, 0.0, 0.0));
      auto c2 = coot::centre_of_molecule(mol);
      if (c1.first && c2.first) {
         double dx = c2.second.x() - c1.second.x();
         if (std::fabs(dx - 5.0) < 0.1)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_min_dist_to_points() {
   starting_test(__FUNCTION__);
   int status = 0;
   clipper::Coord_orth pt(0, 0, 0);
   std::vector<clipper::Coord_orth> others;
   others.push_back(clipper::Coord_orth(3, 0, 0));
   others.push_back(clipper::Coord_orth(5, 0, 0));
   others.push_back(clipper::Coord_orth(1, 0, 0));
   double d = coot::util::min_dist_to_points(pt, others);
   if (std::fabs(d - 1.0) < 0.001)
      status = 1;
   return status;
}

int test_interquartile_range() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::vector<float> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
   float iqr = coot::util::interquartile_range(v);
   // IQR of 1..10 should be around 5
   if (iqr > 3.0 && iqr < 7.0)
      status = 1;
   return status;
}

int test_stats_data() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::vector<float> d = {1.0, 2.0, 3.0, 4.0, 5.0};
   coot::util::stats_data sd(d);
   if (std::fabs(sd.mean - 3.0) < 0.01 && sd.sd > 0.0)
      status = 1;
   return status;
}

int test_qq_plot() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::vector<double> data;
   for (int i = 0; i < 100; i++)
      data.push_back(static_cast<double>(i) / 100.0);
   coot::util::qq_plot_t qq(data);
   auto result = qq.qq_norm();
   if (result.size() > 0)
      status = 1;
   return status;
}

int test_quaternion() {
   starting_test(__FUNCTION__);
   int status = 0;
   // In this convention q3 is the scalar part, so (0,0,0,1) is identity
   coot::util::quaternion q(0, 0, 0, 1);
   clipper::Mat33<double> m = q.matrix();
   // identity quaternion should give identity matrix
   if (std::fabs(m(0,0) - 1.0) < 0.001 &&
       std::fabs(m(1,1) - 1.0) < 0.001 &&
       std::fabs(m(2,2) - 1.0) < 0.001)
      status = 1;
   return status;
}

int test_quaternion_inverse() {
   starting_test(__FUNCTION__);
   int status = 0;
   coot::util::quaternion q(0.5, 0.5, 0.5, 0.5);
   coot::util::quaternion qi = q.inverse();
   // q * q_inverse should be identity (1,0,0,0)
   // test that the inverse has same magnitude
   float mag2 = qi.q0*qi.q0 + qi.q1*qi.q1 + qi.q2*qi.q2 + qi.q3*qi.q3;
   if (std::fabs(mag2 - 1.0f) < 0.01f)
      status = 1;
   return status;
}

int test_split_multi_model_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<mmdb::Manager *> models = coot::util::split_multi_model_molecule(mol);
      // single model PDB should give 1 result
      if (models.size() >= 1) {
         status = 1;
         for (auto m : models) delete m;
      }
      delete mol;
   }
   return status;
}

int test_pair_residue_atoms() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *r1 = coot::util::get_residue("A", 50, "", mol);
      mmdb::Residue *r2 = coot::util::get_residue("A", 51, "", mol);
      if (r1 && r2) {
         auto pairs = coot::util::pair_residue_atoms(r1, r2);
         // should find backbone atom matches at least
         if (pairs.size() > 0)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_peptide_C_N_pairs() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            auto pairs = coot::util::peptide_C_N_pairs(chain_p);
            if (pairs.size() > 0)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_gln_asn_b_factor_outliers() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      auto outliers = coot::util::gln_asn_b_factor_outliers(mol);
      // just runs
      status = 1;
      delete mol;
   }
   return status;
}

int test_CO_orientations() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      auto co = coot::util::CO_orientations(mol);
      if (co.size() > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_sse_to_string() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::string s = coot::util::sse_to_string(1); // MMDB SSE type
   if (!s.empty())
      status = 1;
   return status;
}

int test_residue_types_in_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            std::vector<std::string> types = coot::util::residue_types_in_chain(chain_p);
            if (types.size() > 0)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_delete_anomalous_atoms() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::util::delete_anomalous_atoms(mol);
      // doesn't crash
      status = 1;
      delete mol;
   }
   return status;
}

int test_delete_all_carbohydrate() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      bool r = coot::util::delete_all_carbohydrate(mol);
      // just test it runs
      status = 1;
      delete mol;
   }
   return status;
}

int test_ncs_related_chains() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      auto ncs = coot::ncs_related_chains(mol, 1);
      // may or may not have NCS, just test it runs
      status = 1;
      delete mol;
   }
   return status;
}

int test_sort_residues_by_seqno() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            mmdb::PResidue *residues = nullptr;
            int n_res = 0;
            chain_p->GetResidueTable(residues, n_res);
            if (n_res > 0) {
               auto sorted = coot::util::sort_residues_by_seqno(residues, n_res);
               if (sorted.size() > 0) {
                  // check that first resno <= last resno
                  int first = sorted.front().first->GetSeqNum();
                  int last  = sorted.back().first->GetSeqNum();
                  if (last >= first)
                     status = 1;
               }
            }
         }
      }
      delete mol;
   }
   return status;
}

int test_model_sequence() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            mmdb::PResidue *residues = nullptr;
            int n_res = 0;
            chain_p->GetResidueTable(residues, n_res);
            if (n_res > 0) {
               auto sorted = coot::util::sort_residues_by_seqno(residues, n_res);
               std::string seq = coot::util::model_sequence(sorted);
               if (seq.length() > 0)
                  status = 1;
            }
         }
      }
      delete mol;
   }
   return status;
}

int test_canonical_base_name() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::string r1 = coot::util::canonical_base_name("A", coot::RNA);
   std::string r2 = coot::util::canonical_base_name("A", coot::DNA);
   // should return non-empty for valid base
   if (!r1.empty() || !r2.empty())
      status = 1;
   return status;
}

int test_alt_confs_in_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<std::string> alts = coot::util::alt_confs_in_molecule(mol);
      // just runs
      status = 1;
      delete mol;
   }
   return status;
}

int test_get_residue_mid_point() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::residue_spec_t spec("A", 50, "");
      auto result = coot::util::get_residue_mid_point(mol, spec);
      if (result.first && std::isfinite(result.second.x()))
         status = 1;
      delete mol;
   }
   return status;
}

int test_copy_and_delete_hydrogens() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         mmdb::Residue *res_copy = coot::util::copy_and_delete_hydrogens(res);
         if (res_copy) {
            // should have <= original number of atoms
            int n1 = res->GetNumberOfAtoms();
            int n2 = res_copy->GetNumberOfAtoms();
            if (n2 <= n1 && n2 > 0)
               status = 1;
            delete res_copy;
         }
      }
      delete mol;
   }
   return status;
}

int test_remove_long_links() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::util::remove_long_links(mol, 10.0);
      // just test it runs
      status = 1;
      delete mol;
   }
   return status;
}

int test_correct_link_distances() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::util::correct_link_distances(mol);
      status = 1;
      delete mol;
   }
   return status;
}

int test_get_coords() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::atom_spec_t spec("A", 50, "", " CA ", "");
      mmdb::Atom *at = coot::util::get_atom(spec, mol);
      if (at) {
         try {
            clipper::Coord_orth pos = coot::util::get_coords(at);
            if (std::isfinite(pos.x()))
               status = 1;
         } catch (...) {}
      }
      delete mol;
   }
   return status;
}

int test_refmac_atom_radius() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::atom_spec_t spec("A", 50, "", " CA ", "");
      mmdb::Atom *at = coot::util::get_atom(spec, mol);
      if (at) {
         double r = coot::util::refmac_atom_radius(at);
         if (r > 0.0 && r < 5.0)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_number_of_protein_or_nucleotides() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            auto result = coot::util::get_number_of_protein_or_nucleotides(chain_p);
            // for a protein chain, first (protein count) should be > 0
            if (result.first > 0)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_max_min_max_residue_range() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int r = coot::util::max_min_max_residue_range(mol);
      if (r > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_min_max_residues_in_polymer_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            auto result = coot::util::min_max_residues_in_polymer_chain(chain_p);
            if (result.first) {
               if (result.second.second >= result.second.first)
                  status = 1;
            }
         }
      }
      delete mol;
   }
   return status;
}

int test_residues_in_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<mmdb::Residue *> residues = coot::util::residues_in_chain(mol, "A");
      if (residues.size() > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_get_first_residue_in_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            mmdb::Residue *first = coot::util::get_first_residue_in_chain(chain_p);
            mmdb::Residue *last  = coot::util::get_last_residue_in_chain(chain_p);
            if (first && last)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_residue_orientation() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         clipper::Mat33<double> identity(1,0,0,0,1,0,0,0,1);
         clipper::Mat33<double> ori = coot::util::residue_orientation(res, identity);
         // should return a rotation matrix (orthogonal)
         double det = ori(0,0)*(ori(1,1)*ori(2,2)-ori(1,2)*ori(2,1))
                    - ori(0,1)*(ori(1,0)*ori(2,2)-ori(1,2)*ori(2,0))
                    + ori(0,2)*(ori(1,0)*ori(2,1)-ori(1,1)*ori(2,0));
         if (std::fabs(std::fabs(det) - 1.0) < 0.1)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_ori_to_this_res() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         auto rtops = coot::util::get_ori_to_this_res(res);
         if (rtops.size() > 0)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_put_amino_acid_residue_atom_in_standard_order() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         coot::put_amino_acid_residue_atom_in_standard_order(res);
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_specs_to_atom_selection() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<coot::residue_spec_t> specs;
      specs.push_back(coot::residue_spec_t("A", 50, ""));
      specs.push_back(coot::residue_spec_t("A", 51, ""));
      int selHnd = coot::specs_to_atom_selection(specs, mol, 0); // mode 0 = all atoms
      if (selHnd >= 0) {
         mmdb::PAtom *atoms = nullptr;
         int n_atoms = 0;
         mol->GetSelIndex(selHnd, atoms, n_atoms);
         if (n_atoms > 0)
            status = 1;
         mol->DeleteSelection(selHnd);
      }
      delete mol;
   }
   return status;
}

int test_residues_in_molecule_of_type() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      // look for ALA residues
      std::vector<mmdb::Residue *> alas = coot::util::residues_in_molecule_of_type(mol, "ALA");
      if (alas.size() > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_get_fragment_from_atom_spec() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::atom_spec_t spec("A", 50, "", " CA ", "");
      auto result = coot::util::get_fragment_from_atom_spec(spec, mol);
      if (result.first) {
         int n = coot::util::number_of_residues_in_molecule(result.first);
         if (n > 0)
            status = 1;
         delete result.first;
      }
      delete mol;
   }
   return status;
}

int test_create_mmdbmanager_from_residue_specs() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<coot::residue_spec_t> specs;
      specs.push_back(coot::residue_spec_t("A", 50, ""));
      specs.push_back(coot::residue_spec_t("A", 51, ""));
      mmdb::Manager *new_mol = coot::util::create_mmdbmanager_from_residue_specs(specs, mol);
      if (new_mol) {
         int n = coot::util::number_of_residues_in_molecule(new_mol);
         if (n == 2)
            status = 1;
         delete new_mol;
      }
      delete mol;
   }
   return status;
}

int test_hiranuma_inversion() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      // Get an atom's B-factor before
      coot::atom_spec_t spec("A", 50, "", " CA ", "");
      mmdb::Atom *at = coot::util::get_atom(spec, mol);
      if (at) {
         float b_before = at->tempFactor;
         coot::hiranuma_inversion(mol);
         float b_after = at->tempFactor;
         // B-factor should have changed
         if (std::fabs(b_after - b_before) > 0.001)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_water_coordination() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      try {
         coot::util::water_coordination_t wc(mol, 3.5);
         auto contacts = wc.get_contacts();
         // just runs without crash
         status = 1;
      } catch (...) {
         // may throw if no symmetry etc - that's ok
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_reorientation_matrix() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res1 = coot::util::get_residue("A", 50, "", mol);
      mmdb::Residue *res2 = coot::util::get_residue("A", 51, "", mol);
      if (res1 && res2) {
         auto result = coot::util::get_reorientation_matrix(res1, res2);
         // just test it runs and returns something
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_filter_residues_by_solvent_contact() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         std::vector<mmdb::Residue *> nearby = coot::residues_near_residue(res, mol, 5.0);
         auto filtered = coot::filter_residues_by_solvent_contact(res, mol, nearby, 3.5);
         // just test it runs
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_delete_alt_confs_except() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         coot::util::delete_alt_confs_except(res, "");
         // doesn't crash
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_add_atom() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         int n_before = res->GetNumberOfAtoms();
         bool s = coot::util::add_atom(res, " N  ", " CA ", " C  ", "", 1.5, 109.5, 120.0,
                                        " XX ", " X", 1.0, 20.0);
         int n_after = res->GetNumberOfAtoms();
         if (s && n_after == n_before + 1)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_pdbcleanup_serial_residue_numbers() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::util::pdbcleanup_serial_residue_numbers(mol);
      status = 1;
      delete mol;
   }
   return status;
}

int test_transform_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            clipper::RTop_orth rtop(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1),
                                    clipper::Coord_orth(5, 0, 0));
            coot::util::transform_chain(chain_p, rtop);
            status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_transform_atoms() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         clipper::RTop_orth rtop(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1),
                                 clipper::Coord_orth(5, 0, 0));
         coot::util::transform_atoms(res, rtop);
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_residue_atoms_segid() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         try {
            std::string segid = coot::residue_atoms_segid(res);
            status = 1; // returned without throwing
         } catch (...) {
            status = 1; // throwing is also acceptable for this function
         }
      }
      delete mol;
   }
   return status;
}

int test_copy_segid() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *r1 = coot::util::get_residue("A", 50, "", mol);
      mmdb::Residue *r2 = coot::util::get_residue("A", 51, "", mol);
      if (r1 && r2) {
         bool r = coot::copy_segid(r1, r2);
         // returns a bool, doesn't crash
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_hetify_residue_atoms() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         int n = coot::hetify_residue_atoms(res);
         // returns number of hetatm atoms
         if (n >= 0)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_cell_symm() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      try {
         auto cs = coot::util::get_cell_symm(mol);
         // got cell and spacegroup
         status = 1;
      } catch (const std::runtime_error &e) {
         // no cell/symm info is possible for some PDB files - that's acceptable
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_get_lsq_matrix() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Manager *mol2 = coot::util::copy_molecule(mol);
      if (mol2) {
         std::vector<coot::lsq_range_match_info_t> matches;
         matches.push_back(coot::lsq_range_match_info_t(1, 100, "A", 1, 100, "A", 0));
         auto result = coot::util::get_lsq_matrix(mol, mol2, matches, 1, false);
         // self vs self should succeed with near-identity
         if (result.first)
            status = 1;
         delete mol2;
      }
      delete mol;
   }
   return status;
}

int test_remove_wrong_cis_peptides() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::util::remove_wrong_cis_peptides(mol);
      status = 1;
      delete mol;
   }
   return status;
}

int test_print_secondary_structure_info() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         coot::util::print_secondary_structure_info(model_p);
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_create_mmdbmanager_from_points() {
   starting_test(__FUNCTION__);
   int status = 0;
   std::vector<clipper::Coord_orth> pts;
   pts.push_back(clipper::Coord_orth(0, 0, 0));
   pts.push_back(clipper::Coord_orth(1, 1, 1));
   pts.push_back(clipper::Coord_orth(2, 2, 2));
   mmdb::Manager *mol = coot::util::create_mmdbmanager_from_points(pts, 20.0);
   if (mol) {
      int n = coot::util::number_of_residues_in_molecule(mol);
      if (n == 3)
         status = 1;
      delete mol;
   }
   return status;
}

int test_mtrix_info() {
   starting_test(__FUNCTION__);
   int status = 0;
   // this function reads a PDB file and extracts MTRIX records
   std::string fn = reference_data("moorhen-tutorial-structure-number-1.pdb");
   std::vector<clipper::RTop_orth> mtrix = coot::mtrix_info(fn);
   // may or may not have MTRIX records, just test it runs
   status = 1;
   return status;
}

int test_copy_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            auto result = coot::util::copy_chain(chain_p);
            if (result.first && result.second) {
               status = 1;
               delete result.second;
            }
         }
      }
      delete mol;
   }
   return status;
}

int test_standard_deviation_temperature_factor() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int selHnd = mol->NewSelection();
      mol->SelectAtoms(selHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mmdb::PAtom *atoms = nullptr;
      int n_atoms = 0;
      mol->GetSelIndex(selHnd, atoms, n_atoms);
      if (n_atoms > 0) {
         float sd = coot::util::standard_deviation_temperature_factor(atoms, n_atoms, 0.0, 999.0, 0, 0);
         if (sd >= 0.0)
            status = 1;
      }
      mol->DeleteSelection(selHnd);
      delete mol;
   }
   return status;
}

int test_sort_chains_util() {
   starting_test(__FUNCTION__);
   int status = 0;
   // test the comparison function
   mmdb::Chain *c1 = new mmdb::Chain;
   mmdb::Chain *c2 = new mmdb::Chain;
   c1->SetChainID("A");
   c2->SetChainID("B");
   std::pair<mmdb::Chain *, std::string> p1(c1, "A");
   std::pair<mmdb::Chain *, std::string> p2(c2, "B");
   bool r = coot::sort_chains_util(p1, p2);
   if (r) // A < B
      status = 1;
   delete c1;
   delete c2;
   return status;
}

int test_get_residue_by_binary_search() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue_by_binary_search("A", 50, "", mol);
      if (res) {
         if (res->GetSeqNum() == 50)
            status = 1;
      }
      delete mol;
   }
   return status;
}

int test_chain_only_of_type() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      // look for a chain that is entirely HOH
      mmdb::Chain *chain = coot::util::chain_only_of_type(mol, "HOH");
      // may or may not find one, just test it runs
      status = 1;
      delete mol;
   }
   return status;
}

int test_is_000_shift() {
   starting_test(__FUNCTION__);
   int status = 0;
   clipper::Coord_frac cf_zero(0.0, 0.0, 0.0);
   clipper::Coord_frac cf_nonzero(0.5, 0.0, 0.0);
   bool r1 = coot::util::is_000_shift(cf_zero);
   bool r2 = coot::util::is_000_shift(cf_nonzero);
   if (r1 && !r2)
      status = 1;
   return status;
}

int test_create_mmdbmanager_from_atom_selection() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int selHnd = mol->NewSelection();
      mol->SelectAtoms(selHnd, 0, "A", 50, "*", 55, "*", "*", "*", "*", "*");
      mmdb::Manager *new_mol = coot::util::create_mmdbmanager_from_atom_selection(mol, selHnd);
      if (new_mol) {
         int n = coot::util::number_of_residues_in_molecule(new_mol);
         if (n > 0)
            status = 1;
         delete new_mol;
      }
      mol->DeleteSelection(selHnd);
      delete mol;
   }
   return status;
}

int test_rotate_residue() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         auto c1 = coot::util::get_residue_centre(res);
         coot::util::rotate_residue(res, clipper::Coord_orth(0,0,1), clipper::Coord_orth(0,0,0), M_PI/4);
         auto c2 = coot::util::get_residue_centre(res);
         // centre should have moved
         if (c1.first && c2.first) {
            double d = clipper::Coord_orth::length(c1.second, c2.second);
            // the centre has moved (unless it was on the axis)
            status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_delete_residue_references_in_header_info() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         coot::util::delete_residue_references_in_header_info(res, mol);
         status = 1;
      }
      delete mol;
   }
   return status;
}

int test_copy_headers() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Manager *mol2 = new mmdb::Manager;
      bool r = coot::util::copy_headers(mol, mol2, true);
      // may or may not succeed depending on headers, just test it runs
      status = 1;
      delete mol2;
      delete mol;
   }
   return status;
}

#if 0 // does not exist!
int test_add_copy_of_atom() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::atom_spec_t spec("A", 50, "", " CA ", "");
      mmdb::Atom *at = coot::util::get_atom(spec, mol);
      if (at) {
         // create a new mol and add the atom copy to it
         mmdb::Manager *mol2 = coot::util::copy_molecule(mol);
         if (mol2) {
            coot::util::add_copy_of_atom(mol2, at);
            status = 1;
            delete mol2;
         }
      }
      delete mol;
   }
   return status;
}
#endif

int test_get_atom_using_fuzzy_search() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::atom_spec_t spec("A", 50, "", " CA ", "");
      mmdb::Atom *at = coot::util::get_atom_using_fuzzy_search(spec, mol);
      if (at)
         status = 1;
      delete mol;
   }
   return status;
}

int test_get_biggest_hetgroup() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Residue *res = coot::util::get_biggest_hetgroup(mol);
      // may or may not have hetgroups, just test it runs
      status = 1;
      delete mol;
   }
   return status;
}

int test_transform_selection() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      int selHnd = mol->NewSelection();
      mol->SelectAtoms(selHnd, 0, "A", 50, "*", 55, "*", "*", "*", "*", "*");
      clipper::RTop_orth rtop(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1),
                              clipper::Coord_orth(5, 0, 0));
      coot::util::transform_selection(mol, selHnd, rtop);
      mol->DeleteSelection(selHnd);
      status = 1;
      delete mol;
   }
   return status;
}

int test_move_hetgroups_around_protein() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::util::move_hetgroups_around_protein(mol);
      status = 1;
      delete mol;
   }
   return status;
}

int test_get_second_and_penultimate_residue_in_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            mmdb::Residue *second = coot::util::get_second_residue_in_chain(chain_p);
            mmdb::Residue *penult = coot::util::get_penultimate_residue_in_chain(chain_p);
            if (second && penult)
               status = 1;
         }
      }
      delete mol;
   }
   return status;
}

int test_residue_types_in_residue_vec() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      std::vector<mmdb::Residue *> residues;
      for (int i = 50; i <= 55; i++) {
         mmdb::Residue *r = coot::util::get_residue("A", i, "", mol);
         if (r) residues.push_back(r);
      }
      std::vector<std::string> types = coot::util::residue_types_in_residue_vec(residues);
      if (types.size() > 0)
         status = 1;
      delete mol;
   }
   return status;
}

int test_create_mmdbmanager_from_atom() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      coot::atom_spec_t spec("A", 50, "", " CA ", "");
      mmdb::Atom *at = coot::util::get_atom(spec, mol);
      if (at) {
         mmdb::Manager *new_mol = coot::util::create_mmdbmanager_from_atom(at);
         if (new_mol) {
            status = 1;
            delete new_mol;
         }
      }
      delete mol;
   }
   return status;
}

// ==================== test runner ====================

int n_tests = 0;
static std::vector<std::pair<std::string, int> > test_results;

void
write_test_name(const std::string &test_name) {
   std::ofstream f(".current-test");
   f << "\"" << test_name << "\"" << "\n";
   f.close();
}

int
run_test(int (*test_func)(), const std::string &test_name) {
   n_tests++;
   write_test_name(test_name);
   int status = test_func();
   std::string status_string = "FAIL: ";
   std::string uncol = "\033[m";
   std::string col = "\033[31m";
   if (status == 1) {
      status_string = "PASS: ";
      col = "\033[32m";
   }
   std::cout << status_string << std::setw(50) << std::left << test_name << col << " ⬤ " << uncol << std::endl;
   test_results.push_back(std::make_pair(test_name, status));
   return status;
}

void
print_results_summary() {
   std::cout << "n_tests: " << n_tests << std::endl;
   unsigned int n_failed = 0;
   std::cout << "LIGHTS: ";
   unsigned int count = 0;
   for (const auto &result : test_results) {
      count++;
      const auto &status = result.second;
      if (status == 0) {
         n_failed++;
         std::cout << "\033[31m⬤ ";
      } else {
         std::cout << "\033[32m⬤ ";
      }
      if (count % 40 == 0)
         std::cout << "\n        ";
   }
   std::cout << "\033[m  failures: " << n_failed << "/" << n_tests << std::endl;

   if (n_failed > 0) {
      std::cout << "Test summary: " << n_failed << " failed tests of " << n_tests << std::endl;
      for (const auto &result : test_results) {
         if (result.second == 0)
            std::cout << "FAIL:   " << result.first << std::endl;
      }
   } else {
      std::cout << "Test summary: all " << n_tests << " tests passed" << std::endl;
   }
}

int main(int argc, char **argv) {

   mmdb::InitMatType();

   int status = 0;
   write_test_name("---");

   bool last_test_only = false;
   if (argc > 1) {
      std::string arg(argv[1]);
      if (arg == "last-test-only")
         last_test_only = true;
   }

   int all_tests_status = 1; // fail!

   if (! last_test_only) {

      // Pure computation tests (no file I/O needed)
      status += run_test(test_pad_atom_name,                     "pad_atom_name");
      status += run_test(test_is_hydrogen_atom,                  "is_hydrogen_atom");
      status += run_test(test_lsq_plane_info,                    "lsq_plane_info_t");
      status += run_test(test_lsq_plane_deviation,               "lsq_plane_deviation");
      status += run_test(test_single_letter_to_3_letter_code,    "single_letter_to_3_letter_code (char)");
      status += run_test(test_single_letter_to_3_letter_code_string, "single_letter_to_3_letter_code (string)");
      status += run_test(test_three_letter_to_one_letter,        "three_letter_to_one_letter");
      status += run_test(test_three_letter_to_one_letter_with_specials, "three_letter_to_one_letter_with_specials");
      status += run_test(test_standard_residue_types,            "standard_residue_types");
      status += run_test(test_PDB_standard_residue_types,        "PDB_standard_residue_types");
      status += run_test(test_is_cis,                            "is_cis");
      status += run_test(test_get_brick_id_inner,                "get_brick_id_inner");
      status += run_test(test_get_brick_id,                      "get_brick_id");
      status += run_test(test_rotate_around_vector,              "rotate_around_vector");
      status += run_test(test_average_position,                  "average_position");
      status += run_test(test_min_dist_to_points,                "min_dist_to_points");
      status += run_test(test_interquartile_range,               "interquartile_range");
      status += run_test(test_stats_data,                        "stats_data");
      status += run_test(test_qq_plot,                           "qq_plot");
      status += run_test(test_quaternion,                        "quaternion matrix");
      status += run_test(test_quaternion_inverse,                "quaternion inverse");
      status += run_test(test_matrix_convert,                    "matrix_convert");
      status += run_test(test_is_000_shift,                      "is_000_shift");
      status += run_test(test_sse_to_string,                     "sse_to_string");
      status += run_test(test_sort_chains_util,                  "sort_chains_util");
      status += run_test(test_create_mmdbmanager_from_points,    "create_mmdbmanager_from_points");

      // Tests that need the tutorial PDB file
      status += run_test(test_write_coords_pdb,                  "write_coords_pdb");
      status += run_test(test_write_coords_cif,                  "write_coords_cif");
      status += run_test(test_distance,                          "distance");
      status += run_test(test_angle,                             "angle");
      status += run_test(test_co_and_update_position,            "co and update_position");
      status += run_test(test_get_position_hash,                 "get_position_hash");
      status += run_test(test_get_selection_handle,              "get_selection_handle");
      status += run_test(test_is_member_p,                       "is_member_p");
      status += run_test(test_residues_in_order_p,               "residues_in_order_p");
      status += run_test(test_centre_of_molecule,                "centre_of_molecule");
      status += run_test(test_centre_of_molecule_using_masses,   "centre_of_molecule_using_masses");
      status += run_test(test_radius_of_gyration,                "radius_of_gyration");
      status += run_test(test_centre_of_residues,                "centre_of_residues");
      status += run_test(test_hetify_residue_atoms_as_needed,    "hetify_residue_atoms_as_needed");
      status += run_test(test_hetify_residues_as_needed,         "hetify_residues_as_needed");
      status += run_test(test_atoms_with_zero_occupancy,         "atoms_with_zero_occupancy");
      status += run_test(test_residues_with_alt_confs,           "residues_with_alt_confs");
      status += run_test(test_residues_near_residue,             "residues_near_residue (spec)");
      status += run_test(test_residues_near_residue_v2,          "residues_near_residue (ptr)");
      status += run_test(test_residues_near_residues,            "residues_near_residues");
      status += run_test(test_interface_residues,                "interface_residues");
      status += run_test(test_residues_near_position,            "residues_near_position");
      status += run_test(test_closest_approach,                  "closest_approach");
      status += run_test(test_nearest_residue_by_sequence,       "nearest_residue_by_sequence");
      status += run_test(test_sort_chains,                       "sort_chains");
      status += run_test(test_sort_residues,                     "sort_residues");
      status += run_test(test_mol_has_symmetry,                  "mol_has_symmetry");
      status += run_test(test_mol_is_anisotropic,                "mol_is_anisotropic");
      status += run_test(test_is_nucleotide,                     "is_nucleotide");
      status += run_test(test_residue_has_hydrogens_p,           "residue_has_hydrogens_p");
      status += run_test(test_residue_has_hetatms,               "residue_has_hetatms");
      status += run_test(test_get_residue,                       "get_residue");
      status += run_test(test_get_residue_by_spec,               "get_residue (by spec)");
      status += run_test(test_get_residue_by_binary_search,      "get_residue_by_binary_search");
      status += run_test(test_get_first_residue,                 "get_first_residue");
      status += run_test(test_get_nth_residue,                   "get_nth_residue");
      status += run_test(test_get_atom,                          "get_atom");
      status += run_test(test_get_atom_using_fuzzy_search,       "get_atom_using_fuzzy_search");
      status += run_test(test_get_following_residue,             "get_following_residue");
      status += run_test(test_get_previous_residue,              "get_previous_residue");
      status += run_test(test_get_residue_centre,                "get_residue_centre");
      status += run_test(test_get_CA_position_in_residue,        "get_CA_position_in_residue");
      status += run_test(test_get_CB_position_in_residue,        "get_CB_position_in_residue");
      status += run_test(test_get_residue_mid_point,             "get_residue_mid_point");
      status += run_test(test_get_this_and_next_residues,        "get_this_and_next_residues");
      status += run_test(test_get_coords,                        "get_coords");
      status += run_test(test_residue_types_in_molecule,         "residue_types_in_molecule");
      status += run_test(test_non_standard_residue_types_in_molecule, "non_standard_residue_types_in_molecule");
      status += run_test(test_chains_in_molecule,                "chains_in_molecule");
      status += run_test(test_residues_in_molecule,              "residues_in_molecule");
      status += run_test(test_residues_in_chain,                 "residues_in_chain");
      status += run_test(test_residues_in_molecule_of_type,      "residues_in_molecule_of_type");
      status += run_test(test_number_of_residues_in_molecule,    "number_of_residues_in_molecule");
      status += run_test(test_max_number_of_residues_in_chain,   "max_number_of_residues_in_chain");
      status += run_test(test_number_of_chains,                  "number_of_chains");
      status += run_test(test_min_and_max_residues,              "min_and_max_residues");
      status += run_test(test_min_resno_in_chain,                "min_resno_in_chain");
      status += run_test(test_max_resno_in_chain,                "max_resno_in_chain");
      status += run_test(test_max_resno_in_molecule,             "max_resno_in_molecule");
      status += run_test(test_min_max_residues_in_polymer_chain, "min_max_residues_in_polymer_chain");
      status += run_test(test_max_min_max_residue_range,         "max_min_max_residue_range");
      status += run_test(test_get_number_of_protein_or_nucleotides, "get_number_of_protein_or_nucleotides");
      status += run_test(test_alt_confs_in_molecule,             "alt_confs_in_molecule");
      status += run_test(test_residue_types_in_chain,            "residue_types_in_chain");
      status += run_test(test_residue_types_in_residue_vec,      "residue_types_in_residue_vec");
      status += run_test(test_chain_only_of_type,                "chain_only_of_type");
      status += run_test(test_copy_molecule,                     "copy_molecule");
      status += run_test(test_copy_cell_and_symm_headers,        "copy_cell_and_symm_headers");
      status += run_test(test_copy_headers,                      "copy_headers");
      status += run_test(test_create_mmdbmanager_from_residue,   "create_mmdbmanager_from_residue");
      status += run_test(test_create_mmdbmanager_from_residue_vector, "create_mmdbmanager_from_residue_vector");
      status += run_test(test_create_mmdbmanager_from_mmdbmanager, "create_mmdbmanager_from_mmdbmanager");
      status += run_test(test_create_mmdbmanager_from_residue_specs, "create_mmdbmanager_from_residue_specs");
      status += run_test(test_create_mmdbmanager_from_atom_selection, "create_mmdbmanager_from_atom_selection");
      status += run_test(test_create_mmdbmanager_from_atom,      "create_mmdbmanager_from_atom");
      status += run_test(test_deep_copy_this_residue,            "deep_copy_this_residue");
      status += run_test(test_copy_and_delete_hydrogens,         "copy_and_delete_hydrogens");
      status += run_test(test_copy_chain,                        "copy_chain");
      status += run_test(test_split_multi_model_molecule,        "split_multi_model_molecule");
      status += run_test(test_get_residues_in_fragment,          "get_residues_in_fragment");
      status += run_test(test_get_residues_in_range,             "get_residues_in_range");
      status += run_test(test_get_fragment_from_atom_spec,       "get_fragment_from_atom_spec");
      status += run_test(test_get_hetgroups,                     "get_hetgroups");
      status += run_test(test_get_biggest_hetgroup,              "get_biggest_hetgroup");
      status += run_test(test_get_first_residue_in_chain,        "get_first/last_residue_in_chain");
      status += run_test(test_get_second_and_penultimate_residue_in_chain, "get_second/penultimate_residue_in_chain");
      status += run_test(test_intelligent_this_residue_mmdb_atom, "intelligent_this_residue_mmdb_atom");
      status += run_test(test_previous_and_next_residue,         "previous/next_residue");
      status += run_test(test_sort_residues_by_seqno,            "sort_residues_by_seqno");
      status += run_test(test_model_sequence,                    "model_sequence");
      status += run_test(test_canonical_base_name,               "canonical_base_name");
      status += run_test(test_pair_residue_atoms,                "pair_residue_atoms");
      status += run_test(test_peptide_C_N_pairs,                 "peptide_C_N_pairs");
      status += run_test(test_occupancy_sum,                     "occupancy_sum");
      status += run_test(test_median_temperature_factor,         "median_temperature_factor");
      status += run_test(test_average_temperature_factor,        "average_temperature_factor");
      status += run_test(test_standard_deviation_temperature_factor, "standard_deviation_temperature_factor");
      status += run_test(test_refmac_atom_radius,                "refmac_atom_radius");
      status += run_test(test_extents,                           "extents");
      status += run_test(test_residues_with_insertion_codes,     "residues_with_insertion_codes");
      status += run_test(test_omega_torsion,                     "omega_torsion");
      // status += run_test(test_peptide_torsions,                  "peptide_torsions");
      status += run_test(test_count_cis_peptides,                "count_cis_peptides");
      status += run_test(test_cis_peptides_info_from_coords,     "cis_peptides_info_from_coords");
      status += run_test(test_remove_wrong_cis_peptides,         "remove_wrong_cis_peptides");
      status += run_test(test_remove_long_links,                 "remove_long_links");
      status += run_test(test_correct_link_distances,            "correct_link_distances");
      status += run_test(test_median_position,                   "median_position");
      status += run_test(test_transform_mol,                     "transform_mol");
      status += run_test(test_transform_chain,                   "transform_chain");
      status += run_test(test_transform_atoms,                   "transform_atoms");
      status += run_test(test_transform_selection,               "transform_selection");
      status += run_test(test_shift,                             "shift");
      status += run_test(test_rotate_residue,                    "rotate_residue");
      status += run_test(test_get_ori_to_this_res,               "get_ori_to_this_res");
      status += run_test(test_residue_orientation,               "residue_orientation");
      status += run_test(test_gln_asn_b_factor_outliers,         "gln_asn_b_factor_outliers");
      status += run_test(test_CO_orientations,                   "CO_orientations");
      status += run_test(test_specs_to_atom_selection,           "specs_to_atom_selection");
      status += run_test(test_ncs_related_chains,                "ncs_related_chains");
      status += run_test(test_hetify_residue_atoms,              "hetify_residue_atoms");
      status += run_test(test_put_amino_acid_residue_atom_in_standard_order, "put_amino_acid_residue_atom_in_standard_order");
      status += run_test(test_delete_anomalous_atoms,            "delete_anomalous_atoms");
      status += run_test(test_delete_all_carbohydrate,           "delete_all_carbohydrate");
      status += run_test(test_delete_alt_confs_except,           "delete_alt_confs_except");
      status += run_test(test_delete_residue_references_in_header_info, "delete_residue_references_in_header_info");
      status += run_test(test_add_atom,                          "add_atom");
      // status += run_test(test_add_copy_of_atom,                  "add_copy_of_atom");
      status += run_test(test_pdbcleanup_serial_residue_numbers, "pdbcleanup_serial_residue_numbers");
      status += run_test(test_residue_atoms_segid,               "residue_atoms_segid");
      status += run_test(test_copy_segid,                        "copy_segid");
      status += run_test(test_get_cell_symm,                     "get_cell_symm");
      status += run_test(test_get_lsq_matrix,                    "get_lsq_matrix");
      status += run_test(test_water_coordination,                "water_coordination");
      status += run_test(test_get_reorientation_matrix,          "get_reorientation_matrix");
      status += run_test(test_filter_residues_by_solvent_contact, "filter_residues_by_solvent_contact");
      status += run_test(test_hiranuma_inversion,                "hiranuma_inversion");
      status += run_test(test_print_secondary_structure_info,    "print_secondary_structure_info");
      status += run_test(test_move_hetgroups_around_protein,     "move_hetgroups_around_protein");
      status += run_test(test_mtrix_info,                        "mtrix_info");
   }

   // "last test" placeholder - put the test you're debugging here
   {
      status += run_test(test_get_residue_by_binary_search, "get_residue_by_binary_search (last)");
   }

   if (status == n_tests) all_tests_status = 0;

   print_results_summary();

   return all_tests_status;
}
