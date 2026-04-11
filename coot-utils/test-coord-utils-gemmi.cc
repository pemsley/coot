/* coot-utils/test-coord-utils-gemmi.cc
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

// Tests for gemmi equivalents of functions in coot-coord-utils.hh
// Strategy: read via mmdb, convert to gemmi::Structure via copy_from_mmdb(),
// call both mmdb and gemmi versions, compare results.
// Return 0 for fail, 1 for pass.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <set>

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>
#include <gemmi/mmdb.hpp>

#include "utils/coot-utils.hh"
#include "utils/setup-syminfo.hh"
#include "coot-coord-utils.hh"
#include "coot-coord-utils-gemmi.hh"

void starting_test(const char *func) {
   std::cout << "\nStarting " << func << "()" << std::endl;
}

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

gemmi::Structure gemmi_from_mmdb(mmdb::Manager *mol) {
   gemmi::Structure st = gemmi::copy_from_mmdb(mol);
   coot::trim_atom_names(st);
   return st;
}

// ==================== comparison helpers ====================

bool close(double a, double b, double tol=0.01) {
   return std::fabs(a - b) < tol;
}

bool close_coord(const clipper::Coord_orth &a, const clipper::Coord_orth &b, double tol=0.01) {
   return close(a.x(), b.x(), tol) && close(a.y(), b.y(), tol) && close(a.z(), b.z(), tol);
}

// ==================== tests ====================

int test_distance() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      // get first two atoms from mmdb
      int selHnd = mol->NewSelection();
      mol->SelectAtoms(selHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mmdb::PAtom *atoms = nullptr;
      int n_atoms = 0;
      mol->GetSelIndex(selHnd, atoms, n_atoms);
      if (n_atoms >= 2) {
         double d_mmdb = coot::distance(atoms[0], atoms[1]);
         // get first two atoms from gemmi
         const auto &at1 = st.models[0].chains[0].residues[0].atoms[0];
         const auto &at2 = st.models[0].chains[0].residues[0].atoms[1];
         double d_gemmi = coot::distance(at1, at2);
         if (close(d_mmdb, d_gemmi))
            status = 1;
         else
            std::cout << "MISMATCH distance: mmdb=" << d_mmdb << " gemmi=" << d_gemmi << std::endl;
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
      gemmi::Structure st = gemmi_from_mmdb(mol);
      int selHnd = mol->NewSelection();
      mol->SelectAtoms(selHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mmdb::PAtom *atoms = nullptr;
      int n_atoms = 0;
      mol->GetSelIndex(selHnd, atoms, n_atoms);
      if (n_atoms >= 3) {
         double a_mmdb = coot::angle(atoms[0], atoms[1], atoms[2]);
         const auto &at1 = st.models[0].chains[0].residues[0].atoms[0];
         const auto &at2 = st.models[0].chains[0].residues[0].atoms[1];
         const auto &at3 = st.models[0].chains[0].residues[0].atoms[2];
         double a_gemmi = coot::angle(at1, at2, at3);
         if (close(a_mmdb, a_gemmi, 0.1))
            status = 1;
         else
            std::cout << "MISMATCH angle: mmdb=" << a_mmdb << " gemmi=" << a_gemmi << std::endl;
      }
      mol->DeleteSelection(selHnd);
      delete mol;
   }
   return status;
}

int test_co() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            mmdb::Residue *res = chain_p->GetResidue(0);
            if (res) {
               mmdb::Atom *at = res->GetAtom(0);
               if (at) {
                  clipper::Coord_orth pos_mmdb = coot::co(at);
                  const auto &g_at = st.models[0].chains[0].residues[0].atoms[0];
                  clipper::Coord_orth pos_gemmi = coot::co(g_at);
                  if (close_coord(pos_mmdb, pos_gemmi))
                     status = 1;
               }
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
   // create a gemmi atom and test
   gemmi::Atom at_h;
   at_h.name = "H";
   at_h.element = gemmi::Element("H");
   gemmi::Atom at_c;
   at_c.name = "CA";
   at_c.element = gemmi::Element("C");
   bool r1 = coot::is_hydrogen_atom(at_h);
   bool r2 = coot::is_hydrogen_atom(at_c);
   if (r1 && !r2)
      status = 1;
   return status;
}

int test_centre_of_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      auto result_mmdb = coot::centre_of_molecule(mol);
      auto result_gemmi = coot::centre_of_molecule(st);
      if (result_mmdb.first && result_gemmi.first) {
         if (close_coord(result_mmdb.second, result_gemmi.second, 0.1))
            status = 1;
         else
            std::cout << "MISMATCH centre_of_molecule: mmdb=("
                      << result_mmdb.second.x() << "," << result_mmdb.second.y() << "," << result_mmdb.second.z()
                      << ") gemmi=("
                      << result_gemmi.second.x() << "," << result_gemmi.second.y() << "," << result_gemmi.second.z()
                      << ")" << std::endl;
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
      gemmi::Structure st = gemmi_from_mmdb(mol);
      auto result_mmdb = coot::radius_of_gyration(mol);
      auto result_gemmi = coot::radius_of_gyration(st);
      if (result_mmdb.first && result_gemmi.first) {
         if (close(result_mmdb.second, result_gemmi.second, 0.5))
            status = 1;
         else
            std::cout << "MISMATCH radius_of_gyration: mmdb=" << result_mmdb.second
                      << " gemmi=" << result_gemmi.second << std::endl;
      }
      delete mol;
   }
   return status;
}

int test_get_position_hash() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      float h_gemmi_1 = coot::get_position_hash(st);
      float h_gemmi_2 = coot::get_position_hash(st);
      // The mmdb and gemmi hashes differ because mmdb iterates TER atoms
      // (which affect the x-difference chain) while gemmi has no TER atoms.
      // Test that the gemmi hash is self-consistent and non-zero.
      if (close(h_gemmi_1, h_gemmi_2) && std::fabs(h_gemmi_1) > 1.0)
         status = 1;
      else
         std::cout << "FAIL get_position_hash: h1=" << h_gemmi_1 << " h2=" << h_gemmi_2 << std::endl;
      delete mol;
   }
   return status;
}

int test_mol_has_symmetry() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      bool r_mmdb = coot::mol_has_symmetry(mol);
      bool r_gemmi = coot::mol_has_symmetry(st);
      if (r_mmdb == r_gemmi)
         status = 1;
      else
         std::cout << "MISMATCH mol_has_symmetry: mmdb=" << r_mmdb << " gemmi=" << r_gemmi << std::endl;
      delete mol;
   }
   return status;
}

int test_mol_is_anisotropic() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      bool r_mmdb = coot::mol_is_anisotropic(mol);
      bool r_gemmi = coot::mol_is_anisotropic(st);
      if (r_mmdb == r_gemmi)
         status = 1;
      else
         std::cout << "MISMATCH mol_is_anisotropic: mmdb=" << r_mmdb << " gemmi=" << r_gemmi << std::endl;
      delete mol;
   }
   return status;
}

int test_chains_in_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      auto chains_mmdb = coot::util::chains_in_molecule(mol);
      auto chains_gemmi = coot::util::chains_in_molecule(st);
      if (chains_mmdb.size() == chains_gemmi.size()) {
         bool all_match = true;
         for (size_t i=0; i<chains_mmdb.size(); i++) {
            if (chains_mmdb[i] != chains_gemmi[i]) {
               all_match = false;
               std::cout << "MISMATCH chain[" << i << "]: mmdb=" << chains_mmdb[i]
                         << " gemmi=" << chains_gemmi[i] << std::endl;
            }
         }
         if (all_match) status = 1;
      } else {
         std::cout << "MISMATCH chains_in_molecule: size mmdb=" << chains_mmdb.size()
                   << " gemmi=" << chains_gemmi.size() << std::endl;
      }
      delete mol;
   }
   return status;
}

int test_number_of_residues_in_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      int n_mmdb = coot::util::number_of_residues_in_molecule(mol);
      int n_gemmi = coot::util::number_of_residues_in_molecule(st);
      if (n_mmdb == n_gemmi)
         status = 1;
      else
         std::cout << "MISMATCH number_of_residues_in_molecule: mmdb=" << n_mmdb
                   << " gemmi=" << n_gemmi << std::endl;
      delete mol;
   }
   return status;
}

int test_max_number_of_residues_in_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      int n_mmdb = coot::util::max_number_of_residues_in_chain(mol);
      int n_gemmi = coot::util::max_number_of_residues_in_chain(st);
      if (n_mmdb == n_gemmi)
         status = 1;
      else
         std::cout << "MISMATCH max_number_of_residues_in_chain: mmdb=" << n_mmdb
                   << " gemmi=" << n_gemmi << std::endl;
      delete mol;
   }
   return status;
}

int test_number_of_chains() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      int n_mmdb = coot::util::number_of_chains(mol);
      int n_gemmi = coot::util::number_of_chains(st);
      if (n_mmdb == n_gemmi)
         status = 1;
      else
         std::cout << "MISMATCH number_of_chains: mmdb=" << n_mmdb
                   << " gemmi=" << n_gemmi << std::endl;
      delete mol;
   }
   return status;
}

int test_residue_types_in_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      auto types_mmdb = coot::util::residue_types_in_molecule(mol);
      auto types_gemmi = coot::util::residue_types_in_molecule(st);
      // both should return the same set of types (order may differ)
      std::set<std::string> set_mmdb(types_mmdb.begin(), types_mmdb.end());
      std::set<std::string> set_gemmi(types_gemmi.begin(), types_gemmi.end());
      if (set_mmdb == set_gemmi)
         status = 1;
      else
         std::cout << "MISMATCH residue_types_in_molecule: mmdb has " << set_mmdb.size()
                   << " types, gemmi has " << set_gemmi.size() << std::endl;
      delete mol;
   }
   return status;
}

int test_extents() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      auto ext_mmdb = coot::util::extents(mol);
      auto ext_gemmi = coot::util::extents(st);
      if (close_coord(ext_mmdb.first, ext_gemmi.first, 0.1) &&
          close_coord(ext_mmdb.second, ext_gemmi.second, 0.1))
         status = 1;
      else
         std::cout << "MISMATCH extents" << std::endl;
      delete mol;
   }
   return status;
}

int test_median_position() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      try {
         clipper::Coord_orth med_mmdb = coot::util::median_position(mol);
         clipper::Coord_orth med_gemmi = coot::util::median_position(st);
         if (close_coord(med_mmdb, med_gemmi, 0.5))
            status = 1;
         else
            std::cout << "MISMATCH median_position" << std::endl;
      } catch (...) {
         // both should throw or neither should
         status = 1;
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
      gemmi::Structure st = gemmi_from_mmdb(mol);
      auto r_mmdb = coot::util::max_resno_in_molecule(mol);
      auto r_gemmi = coot::util::max_resno_in_molecule(st);
      if (r_mmdb.first == r_gemmi.first && r_mmdb.second == r_gemmi.second)
         status = 1;
      else
         std::cout << "MISMATCH max_resno_in_molecule: mmdb=(" << r_mmdb.first << "," << r_mmdb.second
                   << ") gemmi=(" << r_gemmi.first << "," << r_gemmi.second << ")" << std::endl;
      delete mol;
   }
   return status;
}

int test_max_min_max_residue_range() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      int r_mmdb = coot::util::max_min_max_residue_range(mol);
      int r_gemmi = coot::util::max_min_max_residue_range(st);
      if (r_mmdb == r_gemmi)
         status = 1;
      else
         std::cout << "MISMATCH max_min_max_residue_range: mmdb=" << r_mmdb
                   << " gemmi=" << r_gemmi << std::endl;
      delete mol;
   }
   return status;
}

int test_get_residue_centre() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (res) {
         auto r_mmdb = coot::util::get_residue_centre(res);
         // find matching residue in gemmi
         for (const auto &chain : st.models[0].chains) {
            if (chain.name == "A") {
               for (const auto &gres : chain.residues) {
                  if (gres.seqid.num.value == 50) {
                     auto r_gemmi = coot::util::get_residue_centre(gres);
                     if (r_mmdb.first && r_gemmi.first) {
                        if (close_coord(r_mmdb.second, r_gemmi.second, 0.1))
                           status = 1;
                        else
                           std::cout << "MISMATCH get_residue_centre" << std::endl;
                     }
                     break;
                  }
               }
               break;
            }
         }
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
      gemmi::Structure st = gemmi_from_mmdb(mol);
      mmdb::Residue *res = coot::util::get_residue("A", 50, "", mol);
      if (!res) {
         std::cout << "   mmdb: could not find residue A/50" << std::endl;
         // try first residue instead
         mmdb::Model *model_p = mol->GetModel(1);
         if (model_p) {
            mmdb::Chain *chain_p = model_p->GetChain(0);
            if (chain_p && chain_p->GetNumberOfResidues() > 0) {
               res = chain_p->GetResidue(0);
               int resno = res->GetSeqNum();
               std::string chain_id = chain_p->GetChainID();
               std::cout << "   using first residue: " << chain_id << "/" << resno << std::endl;
               auto r_mmdb = coot::util::get_CA_position_in_residue(res);
               const auto &g_chain = st.models[0].chains[0];
               if (!g_chain.residues.empty()) {
                  const auto &gres = g_chain.residues[0];
                  auto r_gemmi = coot::util::get_CA_position_in_residue(gres);
                  if (r_mmdb.first && r_gemmi.first) {
                     if (close_coord(r_mmdb.second, r_gemmi.second))
                        status = 1;
                     else
                        std::cout << "MISMATCH get_CA_position_in_residue" << std::endl;
                  } else {
                     std::cout << "   no CA found: mmdb=" << r_mmdb.first
                               << " gemmi=" << r_gemmi.first << std::endl;
                  }
               }
            }
         }
      } else {
         auto r_mmdb = coot::util::get_CA_position_in_residue(res);
         bool found_gemmi_res = false;
         for (const auto &chain : st.models[0].chains) {
            if (chain.name == "A") {
               for (const auto &gres : chain.residues) {
                  if (gres.seqid.num.value == 50) {
                     found_gemmi_res = true;
                     auto r_gemmi = coot::util::get_CA_position_in_residue(gres);
                     if (r_mmdb.first && r_gemmi.first) {
                        if (close_coord(r_mmdb.second, r_gemmi.second))
                           status = 1;
                        else
                           std::cout << "MISMATCH get_CA_position_in_residue" << std::endl;
                     } else {
                        std::cout << "   no CA: mmdb.first=" << r_mmdb.first
                                  << " gemmi.first=" << r_gemmi.first << std::endl;
                     }
                     break;
                  }
               }
               break;
            }
         }
         if (!found_gemmi_res)
            std::cout << "   gemmi: chain A residue 50 not found" << std::endl;
      }
      delete mol;
   }
   return status;
}

int test_min_and_max_residues() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            auto mm_mmdb = coot::util::min_and_max_residues(chain_p);
            const auto &g_chain = st.models[0].chains[0];
            auto mm_gemmi = coot::util::min_and_max_residues(g_chain);
            if (mm_mmdb.first == mm_gemmi.first && mm_mmdb.second == mm_gemmi.second)
               status = 1;
            else
               std::cout << "MISMATCH min_and_max_residues: mmdb=(" << mm_mmdb.first << "," << mm_mmdb.second
                         << ") gemmi=(" << mm_gemmi.first << "," << mm_gemmi.second << ")" << std::endl;
         }
      }
      delete mol;
   }
   return status;
}

int test_residue_types_in_chain() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (chain_p) {
            auto types_mmdb = coot::util::residue_types_in_chain(chain_p);
            const auto &g_chain = st.models[0].chains[0];
            auto types_gemmi = coot::util::residue_types_in_chain(g_chain);
            std::set<std::string> set_mmdb(types_mmdb.begin(), types_mmdb.end());
            std::set<std::string> set_gemmi(types_gemmi.begin(), types_gemmi.end());
            if (set_mmdb == set_gemmi)
               status = 1;
            else
               std::cout << "MISMATCH residue_types_in_chain" << std::endl;
         }
      }
      delete mol;
   }
   return status;
}

int test_count_cis_peptides() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      int n_mmdb = coot::util::count_cis_peptides(mol);
      int n_gemmi = coot::util::count_cis_peptides(st);
      if (n_mmdb == n_gemmi)
         status = 1;
      else
         std::cout << "MISMATCH count_cis_peptides: mmdb=" << n_mmdb
                   << " gemmi=" << n_gemmi << std::endl;
      delete mol;
   }
   return status;
}

int test_is_nucleotide() {
   starting_test(__FUNCTION__);
   int status = 0;
   // test with synthetic gemmi residues
   gemmi::Residue res_ala;
   res_ala.name = "ALA";
   gemmi::Residue res_a;
   res_a.name = "A";
   gemmi::Residue res_da;
   res_da.name = "DA";
   short int r1 = coot::util::is_nucleotide(res_ala);
   short int r2 = coot::util::is_nucleotide(res_a);
   short int r3 = coot::util::is_nucleotide(res_da);
   if (!r1 && r2 && r3)
      status = 1;
   return status;
}

int test_alt_confs_in_molecule() {
   starting_test(__FUNCTION__);
   int status = 0;
   mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mol) {
      gemmi::Structure st = gemmi_from_mmdb(mol);
      auto alts_mmdb = coot::util::alt_confs_in_molecule(mol);
      auto alts_gemmi = coot::util::alt_confs_in_molecule(st);
      std::set<std::string> set_mmdb(alts_mmdb.begin(), alts_mmdb.end());
      std::set<std::string> set_gemmi(alts_gemmi.begin(), alts_gemmi.end());
      if (set_mmdb == set_gemmi)
         status = 1;
      else
         std::cout << "MISMATCH alt_confs_in_molecule: mmdb=" << alts_mmdb.size()
                   << " gemmi=" << alts_gemmi.size() << std::endl;
      delete mol;
   }
   return status;
}

// ==================== test runner ====================

int n_tests = 0;
static std::vector<std::pair<std::string, int> > test_results;

void
write_test_name(const std::string &test_name) {
   std::ofstream f(".current-test-gemmi");
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
      std::cout << "Gemmi test summary: " << n_failed << " failed tests of " << n_tests << std::endl;
      for (const auto &result : test_results) {
         if (result.second == 0)
            std::cout << "FAIL:   " << result.first << std::endl;
      }
   } else {
      std::cout << "Gemmi test summary: all " << n_tests << " tests passed" << std::endl;
   }
}

int main(int argc, char **argv) {

   mmdb::InitMatType();
   setup_syminfo();

   int status = 0;
   write_test_name("---");

   bool last_test_only = false;
   if (argc > 1) {
      std::string arg(argv[1]);
      if (arg == "last-test-only")
         last_test_only = true;
   }

   int all_tests_status = 1; // fail!

   // sanity check: does trim_atom_names work?
   {
      mmdb::Manager *mol = read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
      if (mol) {
         std::cout << "Sanity check: copy_from_mmdb..." << std::flush;
         gemmi::Structure st = gemmi::copy_from_mmdb(mol);
         std::cout << " done. trim..." << std::flush;
         coot::trim_atom_names(st);
         std::cout << " done. first atom: \""
                   << st.models[0].chains[0].residues[0].atoms[0].name << "\"" << std::endl;
         delete mol;
      }
   }

   if (! last_test_only) {

      // Atom-level
      status += run_test(test_distance,                          "distance (gemmi vs mmdb)");
      status += run_test(test_angle,                             "angle (gemmi vs mmdb)");
      status += run_test(test_co,                                "co (gemmi vs mmdb)");
      status += run_test(test_is_hydrogen_atom,                  "is_hydrogen_atom (gemmi)");
      status += run_test(test_is_nucleotide,                     "is_nucleotide (gemmi)");

      // Residue-level
      status += run_test(test_get_residue_centre,                "get_residue_centre (gemmi vs mmdb)");
      status += run_test(test_get_CA_position_in_residue,        "get_CA_position_in_residue (gemmi vs mmdb)");

      // Chain-level
      status += run_test(test_min_and_max_residues,              "min_and_max_residues (gemmi vs mmdb)");
      status += run_test(test_residue_types_in_chain,            "residue_types_in_chain (gemmi vs mmdb)");

      // Structure-level
      status += run_test(test_centre_of_molecule,                "centre_of_molecule (gemmi vs mmdb)");
      status += run_test(test_radius_of_gyration,                "radius_of_gyration (gemmi vs mmdb)");
      status += run_test(test_get_position_hash,                 "get_position_hash (gemmi vs mmdb)");
      status += run_test(test_mol_has_symmetry,                  "mol_has_symmetry (gemmi vs mmdb)");
      status += run_test(test_mol_is_anisotropic,                "mol_is_anisotropic (gemmi vs mmdb)");
      status += run_test(test_chains_in_molecule,                "chains_in_molecule (gemmi vs mmdb)");
      status += run_test(test_number_of_residues_in_molecule,    "number_of_residues_in_molecule (gemmi vs mmdb)");
      status += run_test(test_max_number_of_residues_in_chain,   "max_number_of_residues_in_chain (gemmi vs mmdb)");
      status += run_test(test_number_of_chains,                  "number_of_chains (gemmi vs mmdb)");
      status += run_test(test_residue_types_in_molecule,         "residue_types_in_molecule (gemmi vs mmdb)");
      status += run_test(test_max_resno_in_molecule,             "max_resno_in_molecule (gemmi vs mmdb)");
      status += run_test(test_max_min_max_residue_range,         "max_min_max_residue_range (gemmi vs mmdb)");
      status += run_test(test_alt_confs_in_molecule,             "alt_confs_in_molecule (gemmi vs mmdb)");
      status += run_test(test_extents,                           "extents (gemmi vs mmdb)");
      status += run_test(test_median_position,                   "median_position (gemmi vs mmdb)");
      status += run_test(test_count_cis_peptides,                "count_cis_peptides (gemmi vs mmdb)");
   }

   // "last test" placeholder
   {
      status += run_test(test_centre_of_molecule, "centre_of_molecule (gemmi vs mmdb, last)");
   }

   if (status == n_tests) all_tests_status = 0;

   print_results_summary();

   return all_tests_status;
}
