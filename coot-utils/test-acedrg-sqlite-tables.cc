/*
 * coot-utils/test-acedrg-sqlite-tables.cc
 *
 * Copyright 2026 by Global Phasing Ltd.
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

// Correctness check for coot::acedrg_sqlite_tables::fill_chemcomp():
// for every monomer in a corpus of CIF files, compare the bond/angle
// values it fills against gemmi::AcedrgTables loaded straight from the
// full ASCII tables (the reference implementation).

#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gemmi/cif.hpp>
#include <gemmi/chemcomp.hpp>
#include <gemmi/acedrg_tables.hpp>

#include "acedrg-sqlite-tables.hh"

namespace {

   // Find the block that actually holds the atom list. AceDRG/Refmac
   // monomer-library CIFs put a "comp_list" header block before the
   // "comp_XXX" block with the atoms; CCD-style files hold everything in
   // a single block. Return the first block with a non-empty
   // _chem_comp_atom.atom_id column.
   gemmi::cif::Block *find_chemcomp_block(gemmi::cif::Document &doc) {
      for (auto &block : doc.blocks)
         if (block.find_values("_chem_comp_atom.atom_id"))
            return &block;
      return nullptr;
   }

   // Build one angle restraint (value/esd = NAN) per neighbour pair around
   // each centre atom, from cc.rt.bonds -- the same neighbour-pair
   // enumeration coot::acedrg_sqlite_tables::fill_chemcomp() uses to seed
   // its prefetch when a molecule arrives with no angle records (CCD-style
   // input); Task 4's dictionary-restraints converter will need the same
   // synthesis step for such input, so this exercises that exact path.
   void synthesize_angles_from_bonds(gemmi::ChemComp &cc) {
      std::map<std::string, size_t> atom_idx = cc.make_atom_index();
      std::vector<std::vector<size_t> > nbs(cc.atoms.size());
      for (const auto &b : cc.rt.bonds) {
         auto it1 = atom_idx.find(b.id1.atom);
         auto it2 = atom_idx.find(b.id2.atom);
         if (it1 == atom_idx.end() || it2 == atom_idx.end()) continue;
         nbs[it1->second].push_back(it2->second);
         nbs[it2->second].push_back(it1->second);
      }
      cc.rt.angles.clear();
      for (size_t centre = 0; centre < cc.atoms.size(); centre++) {
         for (size_t i = 0; i < nbs[centre].size(); i++) {
            for (size_t j = i + 1; j < nbs[centre].size(); j++) {
               gemmi::Restraints::Angle a;
               a.id1 = {1, cc.atoms[nbs[centre][i]].id};
               a.id2 = {1, cc.atoms[centre].id};
               a.id3 = {1, cc.atoms[nbs[centre][j]].id};
               a.value = NAN;
               a.esd   = NAN;
               cc.rt.angles.push_back(a);
            }
         }
      }
   }

   // Copy of cc with all bond values/esds set to NAN. Angle restraints:
   // if cc already has angle records (monomer-library input) keep them
   // and NaN their values/esds; otherwise (CCD-style input, no angle
   // records) synthesize one NaN-valued angle per bonded neighbour pair
   // around each atom, so fill_restraints()/fill_chemcomp() always has
   // real angle topology to fill values into.
   gemmi::ChemComp strip_restraint_values(const gemmi::ChemComp &cc_in) {
      gemmi::ChemComp cc = cc_in;
      for (auto &b : cc.rt.bonds) {
         b.value         = NAN;
         b.esd           = NAN;
         b.value_nucleus = NAN;
         b.esd_nucleus   = NAN;
      }
      if (cc.rt.angles.empty()) {
         synthesize_angles_from_bonds(cc);
      } else {
         for (auto &a : cc.rt.angles) {
            a.value = NAN;
            a.esd   = NAN;
         }
      }
      return cc;
   }

   bool nan_or_close(double a, double b) {
      if (std::isnan(a) && std::isnan(b)) return true;
      return std::fabs(a - b) < 1e-9;
   }

   std::string angle_key(const gemmi::Restraints::Angle &a) {
      std::string lo = a.id1.atom;
      std::string hi = a.id3.atom;
      if (lo > hi) std::swap(lo, hi);
      return a.id2.atom + "|" + lo + "|" + hi;
   }

   // Compare cc_ref (filled by the reference ASCII-table engine) against
   // cc_test (filled by the SQLite-backed engine under test). Returns
   // empty string on success, an error message otherwise.
   std::string compare_chemcomps(const gemmi::ChemComp &cc_ref, const gemmi::ChemComp &cc_test) {

      if (cc_ref.rt.bonds.size() != cc_test.rt.bonds.size())
         return "bond count mismatch: ref " + std::to_string(cc_ref.rt.bonds.size()) +
                " test " + std::to_string(cc_test.rt.bonds.size());

      std::map<std::string, const gemmi::Restraints::Bond *> test_bonds;
      for (const auto &b : cc_test.rt.bonds)
         test_bonds[b.lexicographic_str()] = &b;

      for (const auto &b_ref : cc_ref.rt.bonds) {
         std::string key = b_ref.lexicographic_str();
         auto it = test_bonds.find(key);
         if (it == test_bonds.end())
            return "bond " + key + " present in reference but not in test";
         const gemmi::Restraints::Bond &b_test = *it->second;
         if (! nan_or_close(b_ref.value, b_test.value))
            return "bond " + key + " value mismatch: ref " + std::to_string(b_ref.value) +
                   " test " + std::to_string(b_test.value);
         if (! nan_or_close(b_ref.esd, b_test.esd))
            return "bond " + key + " esd mismatch: ref " + std::to_string(b_ref.esd) +
                   " test " + std::to_string(b_test.esd);
      }

      if (cc_ref.rt.angles.size() != cc_test.rt.angles.size())
         return "angle count mismatch: ref " + std::to_string(cc_ref.rt.angles.size()) +
                " test " + std::to_string(cc_test.rt.angles.size());

      // Every multi-atom monomer should have picked up real angle
      // topology from synthesize_angles_from_bonds() (or from its own
      // pre-existing angle records) -- a monomer with zero angles here
      // (other than the lone-atom ZN) means the angle path silently
      // produced nothing, which is exactly the gap this check exists to
      // catch.
      if (cc_ref.rt.angles.empty() && cc_ref.atoms.size() > 1 && cc_ref.name != "ZN")
         return "expected nonzero angles for multi-atom monomer " + cc_ref.name + ", got 0";

      std::map<std::string, const gemmi::Restraints::Angle *> test_angles;
      for (const auto &a : cc_test.rt.angles)
         test_angles[angle_key(a)] = &a;

      for (const auto &a_ref : cc_ref.rt.angles) {
         std::string key = angle_key(a_ref);
         auto it = test_angles.find(key);
         if (it == test_angles.end())
            return "angle " + key + " present in reference but not in test";
         const gemmi::Restraints::Angle &a_test = *it->second;
         if (! nan_or_close(a_ref.value, a_test.value))
            return "angle " + key + " value mismatch: ref " + std::to_string(a_ref.value) +
                   " test " + std::to_string(a_test.value);
         if (! nan_or_close(a_ref.esd, a_test.esd))
            return "angle " + key + " esd mismatch: ref " + std::to_string(a_ref.esd) +
                   " test " + std::to_string(a_test.esd);
      }

      return std::string();
   }

   int run_corpus(const std::string &full_ascii_dir,
                  const std::string &data_dir,
                  const std::string &corpus_dir) {

      std::cout << "loading reference ASCII tables from " << full_ascii_dir << " ... " << std::flush;
      gemmi::AcedrgTables ref;
      try {
         ref.load_tables(full_ascii_dir);
      }
      catch (const std::exception &e) {
         std::cout << "\nFAIL: could not load reference tables: " << e.what() << std::endl;
         return 1;
      }
      std::cout << "done" << std::endl;

      coot::acedrg_sqlite_tables t;
      if (! t.init(data_dir)) {
         std::cout << "FAIL: could not init acedrg_sqlite_tables from " << data_dir << std::endl;
         return 1;
      }

      std::vector<std::filesystem::path> cif_paths;
      for (const auto &entry : std::filesystem::directory_iterator(corpus_dir)) {
         if (! entry.is_regular_file()) continue;
         if (entry.path().extension() != ".cif") continue;
         std::string name = entry.path().filename().string();
         // "reference_*.cif" files in this corpus are chem_link outputs
         // from a different test (link generation), not monomers.
         if (name.rfind("reference_", 0) == 0) continue;
         cif_paths.push_back(entry.path());
      }
      std::sort(cif_paths.begin(), cif_paths.end());

      bool all_ok = true;
      int n_tested = 0;

      for (const auto &path : cif_paths) {
         std::string name = path.filename().string();
         try {
            gemmi::cif::Document doc = gemmi::cif::read_file(path.string());
            gemmi::cif::Block *block = find_chemcomp_block(doc);
            if (! block) {
               std::cout << "FAIL " << name << ": no block with _chem_comp_atom.atom_id found" << std::endl;
               all_ok = false;
               continue;
            }
            gemmi::ChemComp cc_orig = gemmi::make_chemcomp_from_block(*block);
            std::string comp_id = cc_orig.name;

            gemmi::ChemComp cc_ref  = strip_restraint_values(cc_orig);
            gemmi::ChemComp cc_test = strip_restraint_values(cc_orig);

            ref.fill_restraints(cc_ref);

            bool status = t.fill_chemcomp(cc_test);
            if (! status) {
               std::cout << "FAIL " << comp_id << ": fill_chemcomp() returned false" << std::endl;
               all_ok = false;
               continue;
            }

            ++n_tested;

            std::string err = compare_chemcomps(cc_ref, cc_test);
            if (err.empty()) {
               std::cout << "PASS " << comp_id << " (" << cc_ref.rt.bonds.size() << " bonds, "
                         << cc_ref.rt.angles.size() << " angles)" << std::endl;
            } else {
               std::cout << "FAIL " << comp_id << ": " << err << std::endl;
               all_ok = false;
            }
         }
         catch (const std::exception &e) {
            std::cout << "FAIL " << name << ": exception: " << e.what() << std::endl;
            all_ok = false;
         }
      }

      std::cout << n_tested << " monomer(s) tested" << std::endl;
      return all_ok ? 0 : 1;
   }

} // namespace

int main(int argc, char **argv) {

   std::vector<std::string> args;
   for (int i = 1; i < argc; i++)
      args.push_back(argv[i]);

   if (args.size() == 4 && args[0] == "corpus")
      return run_corpus(args[1], args[2], args[3]);

   std::cout << "Usage: test-acedrg-sqlite-tables corpus <full-ascii-tables-dir> <data-dir> <corpus-dir>"
             << std::endl;
   return 2;
}
