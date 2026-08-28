/*
 * coot-utils/acedrg-sqlite-tables.cc
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

// SQLite backend for the heavy AceDRG bond and angle tables.
//
// coot-make-acedrg-sqlite (the builder half implemented in this file)
// converts a directory of AceDRG ASCII tables into a small directory of
// cheap tables plus a single acedrg.sqlite file holding the heavy
// bond/angle tables. The reader half (init/is_usable/fill_chemcomp/
// make_bond_and_angle_restraints) is implemented in a later task.
//
// This is a direct port of gemmi's src/acedrg_tables_db.cpp (branch
// drg-tables-sqlite) build_acedrg_sqlite() and its helpers -- schema,
// column order, parsing and pragmas kept identical. See
// 2026-08-14-acedrg-sqlite-tables-in-coot-design.md.

#include "acedrg-sqlite-tables.hh"

#include "utils/xdg-base.hh"

#include <filesystem>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <tuple>
#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <algorithm>

#ifdef USE_SQLITE3
#include <sqlite3.h>
#endif

#include <gemmi/acedrg_tables.hpp>
#include <gemmi/chemcomp.hpp>

// ------------------------------------------------------------------------
// impl
// ------------------------------------------------------------------------

class coot::acedrg_sqlite_tables::impl {
public:
   gemmi::AcedrgTables tables;
#ifdef USE_SQLITE3
   sqlite3 *db = nullptr;
#endif
   bool usable = false;
   ~impl() {
#ifdef USE_SQLITE3
      if (db) sqlite3_close(db);
#endif
   }
#ifdef USE_SQLITE3
   // Fetch this molecule's bond/angle rows from acedrg.sqlite into
   // tables' on-demand caches. Direct port of gemmi's (fork branch
   // drg-tables-sqlite) AcedrgTables::prefetch_for_molecule().
   void prefetch_for_molecule(const std::vector<std::tuple<int, int, std::string> > &bond_keys,
                              const std::vector<std::tuple<int, int, int> > &angle_triples);
#endif
};

coot::acedrg_sqlite_tables::acedrg_sqlite_tables() : pimpl(new impl) {}

coot::acedrg_sqlite_tables::~acedrg_sqlite_tables() = default;

std::string
coot::acedrg_sqlite_tables::default_data_dir() {

   xdg_t xdg;
   std::filesystem::path p = xdg.get_cache_home();
   p /= "acedrg-tables";
   return p.string();
}

#ifdef USE_SQLITE3

bool
coot::acedrg_sqlite_tables::init(const std::string &acedrg_data_dir) {

   std::string sqlite_path = acedrg_data_dir + "/acedrg.sqlite";
   if (! std::filesystem::exists(sqlite_path))
      return false;

   try {
      pimpl->tables.load_tables(acedrg_data_dir);
   }
   catch (const std::exception &e) {
      std::cout << "acedrg_sqlite_tables::init(): " << e.what() << std::endl;
      return false;
   }

   if (sqlite3_open_v2(sqlite_path.c_str(), &pimpl->db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) {
      std::cout << "acedrg_sqlite_tables::init(): cannot open " << sqlite_path
                << ": " << (pimpl->db ? sqlite3_errmsg(pimpl->db) : "unknown") << std::endl;
      if (pimpl->db) { sqlite3_close(pimpl->db); pimpl->db = nullptr; }
      return false;
   }
   sqlite3_exec(pimpl->db, "PRAGMA query_only = ON;", nullptr, nullptr, nullptr);

   pimpl->usable = true;
   return true;
}

bool
coot::acedrg_sqlite_tables::init() {

   return init(default_data_dir());
}

bool
coot::acedrg_sqlite_tables::is_usable() const {

   return pimpl->usable;
}

namespace {

   inline std::string sqlite_text(sqlite3_stmt *st, int col) {
      const unsigned char *s = sqlite3_column_text(st, col);
      return s ? std::string(reinterpret_cast<const char *>(s)) : std::string();
   }

   // SQLite-quote a string value so it can be spliced directly into a
   // WHERE-clause VALUES literal (doubling embedded single quotes).
   std::string sql_quote(const std::string &s) {
      std::string out; out.reserve(s.size() + 2);
      out += '\'';
      for (char c : s) { if (c == '\'') out += "''"; else out += c; }
      out += '\'';
      return out;
   }

} // namespace

void
coot::acedrg_sqlite_tables::impl::prefetch_for_molecule(const std::vector<std::tuple<int, int, std::string> > &bond_keys,
                                                         const std::vector<std::tuple<int, int, int> > &angle_triples) {

   if (! db) return;

   // Clear the on-demand caches before refilling for this molecule.
   tables.bond_idx_1d_.clear();
   tables.bond_idx_full_.clear();
   tables.bond_idx_2d_.clear();
   tables.bond_2d_hybr_keys_.clear();
   tables.bond_full_4prefix_keys_.clear();
   tables.angle_idx_1d_.clear();
   tables.angle_idx_2d_.clear();
   tables.angle_idx_3d_.clear();
   tables.angle_idx_4d_.clear();
   tables.angle_idx_5d_.clear();
   tables.angle_idx_6d_.clear();

   // Dedup the keys.
   std::set<std::tuple<int, int, std::string> > uniq_bonds(bond_keys.begin(), bond_keys.end());
   std::set<std::tuple<int, int, int> > uniq_triples(angle_triples.begin(), angle_triples.end());

   // ---- bond rows: (ha1, ha2, hybr_comb) in {keys}. -----------------------
   // Note: no in_ring filter -- fill_bond's Y/N ring fallback needs both
   // values present.
   if (! uniq_bonds.empty()) {
      std::string sql =
         "SELECT ha1, ha2, hybr_comb, in_ring,"
         " a1_nb2, a2_nb2, a1_nb, a2_nb,"
         " a1_type_m, a2_type_m, a1_type_f, a2_type_f,"
         " value, sigma, count, value_1d, sigma_1d, count_1d"
         " FROM bond_entries WHERE (ha1, ha2, hybr_comb) IN (VALUES ";
      bool first = true;
      for (const auto &t : uniq_bonds) {
         if (! first) sql += ',';
         sql += '(';
         sql += std::to_string(std::get<0>(t));
         sql += ',';
         sql += std::to_string(std::get<1>(t));
         sql += ',';
         sql += sql_quote(std::get<2>(t));
         sql += ')';
         first = false;
      }
      sql += ')';

      sqlite3_stmt *st = nullptr;
      if (sqlite3_prepare_v2(db, sql.c_str(), -1, &st, nullptr) != SQLITE_OK)
         throw std::runtime_error(std::string("acedrg-db: prepare bond: ") + sqlite3_errmsg(db));

      std::string key_buf;  key_buf.reserve(512);
      std::string hybr_buf; hybr_buf.reserve(64);
      while (sqlite3_step(st) == SQLITE_ROW) {
         int ha1 = sqlite3_column_int(st, 0);
         int ha2 = sqlite3_column_int(st, 1);
         std::string hybr_comb = sqlite_text(st, 2);
         std::string in_ring   = sqlite_text(st, 3);
         std::string a1_nb2    = sqlite_text(st, 4);
         std::string a2_nb2    = sqlite_text(st, 5);
         std::string a1_nb     = sqlite_text(st, 6);
         std::string a2_nb     = sqlite_text(st, 7);
         std::string a1_type_m = sqlite_text(st, 8);
         std::string a2_type_m = sqlite_text(st, 9);
         std::string a1_type_f = sqlite_text(st, 10);
         std::string a2_type_f = sqlite_text(st, 11);
         double value    = sqlite3_column_double(st, 12);
         double sigma    = sqlite3_column_double(st, 13);
         int    count    = sqlite3_column_int   (st, 14);
         double value_1d = sqlite3_column_double(st, 15);
         double sigma_1d = sqlite3_column_double(st, 16);
         int    count_1d = sqlite3_column_int   (st, 17);
         tables.insert_bond_row(ha1, ha2, hybr_comb, in_ring,
                                a1_nb2, a2_nb2, a1_nb, a2_nb,
                                a1_type_m, a2_type_m, a1_type_f, a2_type_f,
                                gemmi::CodStats(value, sigma, count),
                                gemmi::CodStats(value_1d, sigma_1d, count_1d),
                                key_buf, hybr_buf);
      }
      sqlite3_finalize(st);
   }

   // ---- angle rows: (ha1, ha2, ha3) in {triples}. -------------------------
   if (! uniq_triples.empty()) {
      std::string sql =
         "SELECT ha1, ha2, ha3, value_key,"
         " a1_root, a2_root, a3_root,"
         " a1_nb2, a2_nb2, a3_nb2,"
         " a1_nb, a2_nb, a3_nb,"
         " a1_type, a2_type, a3_type,"
         " v1,s1,c1, v2,s2,c2, v3,s3,c3, v4,s4,c4, v5,s5,c5, v6,s6,c6"
         " FROM angle_entries WHERE (ha1, ha2, ha3) IN (VALUES ";
      bool first = true;
      for (const auto &t : uniq_triples) {
         if (! first) sql += ',';
         sql += '(';
         sql += std::to_string(std::get<0>(t)); sql += ',';
         sql += std::to_string(std::get<1>(t)); sql += ',';
         sql += std::to_string(std::get<2>(t));
         sql += ')';
         first = false;
      }
      sql += ')';

      sqlite3_stmt *st = nullptr;
      if (sqlite3_prepare_v2(db, sql.c_str(), -1, &st, nullptr) != SQLITE_OK)
         throw std::runtime_error(std::string("acedrg-db: prepare angle: ") + sqlite3_errmsg(db));

      std::string key_buf; key_buf.reserve(1024);
      while (sqlite3_step(st) == SQLITE_ROW) {
         int ha1 = sqlite3_column_int(st, 0);
         int ha2 = sqlite3_column_int(st, 1);
         int ha3 = sqlite3_column_int(st, 2);
         std::string value_key = sqlite_text(st, 3);
         std::string a1_root = sqlite_text(st, 4);
         std::string a2_root = sqlite_text(st, 5);
         std::string a3_root = sqlite_text(st, 6);
         std::string a1_nb2  = sqlite_text(st, 7);
         std::string a2_nb2  = sqlite_text(st, 8);
         std::string a3_nb2  = sqlite_text(st, 9);
         std::string a1_nb   = sqlite_text(st, 10);
         std::string a2_nb   = sqlite_text(st, 11);
         std::string a3_nb   = sqlite_text(st, 12);
         std::string a1_type = sqlite_text(st, 13);
         std::string a2_type = sqlite_text(st, 14);
         std::string a3_type = sqlite_text(st, 15);
         double v[6]; double s[6]; int c[6];
         int col = 16;
         for (int i = 0; i < 6; ++i) {
            v[i] = sqlite3_column_double(st, col++);
            s[i] = sqlite3_column_double(st, col++);
            c[i] = sqlite3_column_int   (st, col++);
         }
         tables.insert_angle_row(ha1, ha2, ha3, value_key,
                                 a1_root, a2_root, a3_root,
                                 a1_nb2, a2_nb2, a3_nb2,
                                 a1_nb,  a2_nb,  a3_nb,
                                 a1_type, a2_type, a3_type,
                                 v, s, c, key_buf);
      }
      sqlite3_finalize(st);
   }
}

bool
coot::acedrg_sqlite_tables::fill_chemcomp(gemmi::ChemComp &cc) {

   if (! is_usable()) return false;
   gemmi::AcedrgTables &tables = pimpl->tables;

   try {
      std::vector<gemmi::CodAtomInfo> atom_info = tables.classify_atoms(cc);
      std::map<std::string, size_t> atom_idx = cc.make_atom_index();

      std::vector<std::tuple<int, int, std::string> > bond_keys;
      bond_keys.reserve(cc.rt.bonds.size());
      for (const auto &b : cc.rt.bonds) {
         auto it1 = atom_idx.find(b.id1.atom);
         auto it2 = atom_idx.find(b.id2.atom);
         if (it1 == atom_idx.end() || it2 == atom_idx.end()) continue;
         const gemmi::CodAtomInfo &ai1 = atom_info[it1->second];
         const gemmi::CodAtomInfo &ai2 = atom_info[it2->second];
         int h1 = ai1.hashing_value;
         int h2 = ai2.hashing_value;
         std::string s1 = gemmi::hybridization_to_string(ai1.hybrid);
         std::string s2 = gemmi::hybridization_to_string(ai2.hybrid);
         // schema convention: hashes ascending; hybridization strings sorted
         // lexicographically and joined with '_'
         if (h1 > h2) std::swap(h1, h2);
         if (s1 > s2) std::swap(s1, s2);
         bond_keys.emplace_back(h1, h2, s1 + "_" + s2);
      }

      std::vector<std::tuple<int, int, int> > angle_triples;
      angle_triples.reserve(cc.rt.angles.size());
      for (const auto &a : cc.rt.angles) {
         auto i1 = atom_idx.find(a.id1.atom);
         auto i2 = atom_idx.find(a.id2.atom);
         auto i3 = atom_idx.find(a.id3.atom);
         if (i1 == atom_idx.end() || i2 == atom_idx.end() || i3 == atom_idx.end()) continue;
         int h1 = atom_info[i1->second].hashing_value;
         int h2 = atom_info[i2->second].hashing_value;
         int h3 = atom_info[i3->second].hashing_value;
         // schema convention: ha1 = centre hash; outer hashes sorted
         angle_triples.emplace_back(h2, std::min(h1, h3), std::max(h1, h3));
      }
      if (angle_triples.empty()) {
         // no angle records yet (CCD-style input): the fill pipeline will
         // derive angles from neighbour pairs around each atom - seed the
         // prefetch with those triples
         std::vector<std::vector<size_t> > nbs(cc.atoms.size());
         for (const auto &b : cc.rt.bonds) {
            auto it1 = atom_idx.find(b.id1.atom);
            auto it2 = atom_idx.find(b.id2.atom);
            if (it1 == atom_idx.end() || it2 == atom_idx.end()) continue;
            nbs[it1->second].push_back(it2->second);
            nbs[it2->second].push_back(it1->second);
         }
         for (size_t centre = 0; centre < cc.atoms.size(); centre++) {
            int hc = atom_info[centre].hashing_value;
            for (size_t i = 0; i < nbs[centre].size(); i++) {
               for (size_t j = i + 1; j < nbs[centre].size(); j++) {
                  int ha = atom_info[nbs[centre][i]].hashing_value;
                  int hb = atom_info[nbs[centre][j]].hashing_value;
                  angle_triples.emplace_back(hc, std::min(ha, hb), std::max(ha, hb));
               }
            }
         }
      }

      pimpl->prefetch_for_molecule(bond_keys, angle_triples);
      tables.fill_restraints(cc);
   }
   catch (const std::exception &e) {
      std::cout << "WARNING:: acedrg_sqlite_tables: " << e.what() << std::endl;
      return false;
   }
   return true;
}

#else // !USE_SQLITE3

bool
coot::acedrg_sqlite_tables::init(const std::string &acedrg_data_dir) {

   return false;
}

bool
coot::acedrg_sqlite_tables::init() {

   return false;
}

bool
coot::acedrg_sqlite_tables::is_usable() const {

   return false;
}

bool
coot::acedrg_sqlite_tables::fill_chemcomp(gemmi::ChemComp &cc) {

   return false;
}

#endif // USE_SQLITE3

namespace {

   // Build a gemmi ChemComp (atoms + typed bonds, no values) from a coot
   // dictionary -- the input to AceDRG atom typing. Also synthesizes one
   // angle record (value/esd = NaN) per bonded-neighbour pair around each
   // atom: gemmi's fill_restraints() only fills *existing* angle records,
   // it never derives them from bond adjacency, so the angle topology has
   // to be supplied up front here or make_bond_and_angle_restraints()
   // would always come back with zero angles. Mirrors
   // synthesize_angles_from_bonds() in test-acedrg-sqlite-tables.cc.
   gemmi::ChemComp
   chemcomp_from_dictionary(const coot::dictionary_residue_restraints_t &rest) {

      gemmi::ChemComp cc;
      cc.name = rest.residue_info.comp_id;
      for (const auto &a : rest.atom_info) {
         gemmi::ChemComp::Atom atom;
         atom.id = a.atom_id;
         atom.el = gemmi::Element(a.type_symbol.c_str()).elem;
         if (a.formal_charge.first)
            atom.charge = static_cast<float>(a.formal_charge.second);
         atom.chem_type = a.type_energy;
         cc.atoms.push_back(atom);
      }
      double nan = std::numeric_limits<double>::quiet_NaN();
      for (const auto &b : rest.bond_restraint) {
         gemmi::Restraints::Bond bond;
         bond.id1 = gemmi::Restraints::AtomId{1, b.atom_id_1()};
         bond.id2 = gemmi::Restraints::AtomId{1, b.atom_id_2()};
         bond.type = gemmi::bond_type_from_string(b.type());
         bond.aromatic = gemmi::is_aromatic_or_deloc(bond.type);
         bond.value = nan; bond.esd = nan;
         bond.value_nucleus = nan; bond.esd_nucleus = nan;
         cc.rt.bonds.push_back(bond);
      }

      // synthesize angle topology from bond adjacency (see comment above)
      std::map<std::string, size_t> atom_idx = cc.make_atom_index();
      std::vector<std::vector<size_t> > nbs(cc.atoms.size());
      for (const auto &b : cc.rt.bonds) {
         auto it1 = atom_idx.find(b.id1.atom);
         auto it2 = atom_idx.find(b.id2.atom);
         if (it1 == atom_idx.end() || it2 == atom_idx.end()) continue;
         nbs[it1->second].push_back(it2->second);
         nbs[it2->second].push_back(it1->second);
      }
      for (std::size_t centre = 0; centre < cc.atoms.size(); centre++) {
         for (std::size_t i = 0; i < nbs[centre].size(); i++) {
            for (std::size_t j = i + 1; j < nbs[centre].size(); j++) {
               gemmi::Restraints::Angle ang;
               ang.id1 = gemmi::Restraints::AtomId{1, cc.atoms[nbs[centre][i]].id};
               ang.id2 = gemmi::Restraints::AtomId{1, cc.atoms[centre].id};
               ang.id3 = gemmi::Restraints::AtomId{1, cc.atoms[nbs[centre][j]].id};
               ang.value = nan;
               ang.esd   = nan;
               cc.rt.angles.push_back(ang);
            }
         }
      }

      return cc;
   }

} // namespace

std::pair<bool, coot::dictionary_residue_restraints_t>
coot::acedrg_sqlite_tables::make_bond_and_angle_restraints(const dictionary_residue_restraints_t &restraints_in) {

   dictionary_residue_restraints_t r(restraints_in.residue_info.comp_id, 1);
   if (! is_usable())
      return std::make_pair(false, r);

   gemmi::ChemComp cc = chemcomp_from_dictionary(restraints_in);
   if (! fill_chemcomp(cc))
      return std::make_pair(false, r);

   // conservatively_replace_with() swaps in whole restraint objects, so
   // carry the input's bond-type strings across
   std::map<std::string, std::string> bond_type_for_pair;
   for (const auto &b : restraints_in.bond_restraint)
      bond_type_for_pair[gemmi::Restraints::lexicographic_str(b.atom_id_1(), b.atom_id_2())] = b.type();

   unsigned int n_filled = 0;
   for (const auto &b : cc.rt.bonds) {
      if (std::isnan(b.value) || std::isnan(b.esd)) continue; // no table hit: keep fallback value
      std::string type = "single";
      auto it = bond_type_for_pair.find(gemmi::Restraints::lexicographic_str(b.id1.atom, b.id2.atom));
      if (it != bond_type_for_pair.end()) type = it->second;
      dict_bond_restraint_t br(b.id1.atom, b.id2.atom, type, b.value, b.esd, 0.0, 0.0, false);
      r.bond_restraint.push_back(br);
      n_filled++;
   }
   for (const auto &a : cc.rt.angles) {
      if (std::isnan(a.value) || std::isnan(a.esd)) continue;
      dict_angle_restraint_t ar(a.id1.atom, a.id2.atom, a.id3.atom, a.value, a.esd);
      r.angle_restraint.push_back(ar);
      n_filled++;
   }
   return std::make_pair(n_filled > 0, r);
}

// ------------------------------------------------------------------------
// build()
// ------------------------------------------------------------------------

#ifdef USE_SQLITE3

namespace {

   // Minimal RAII wrappers around sqlite3 and sqlite3_stmt.
   struct SqliteDB {
      sqlite3 *db = nullptr;
      ~SqliteDB() { if (db) sqlite3_close(db); }
      void open(const std::string &path, int flags) {
         if (sqlite3_open_v2(path.c_str(), &db, flags, nullptr) != SQLITE_OK) {
            std::string err = db ? sqlite3_errmsg(db) : "unknown";
            throw std::runtime_error("acedrg-db: cannot open " + path + ": " + err);
         }
      }
      void exec(const char *sql) {
         char *err = nullptr;
         if (sqlite3_exec(db, sql, nullptr, nullptr, &err) != SQLITE_OK) {
            std::string e = err ? err : "unknown";
            sqlite3_free(err);
            throw std::runtime_error("acedrg-db: SQL failed: " + e + " (" + sql + ")");
         }
      }
   };

   struct SqliteStmt {
      sqlite3_stmt *st = nullptr;
      ~SqliteStmt() { if (st) sqlite3_finalize(st); }
      void prepare(sqlite3 *db, const char *sql) {
         if (sqlite3_prepare_v2(db, sql, -1, &st, nullptr) != SQLITE_OK)
            throw std::runtime_error(std::string("acedrg-db: prepare failed: ") +
                                     sqlite3_errmsg(db));
      }
   };

   // Skip whitespace, return pointer to next non-blank.
   inline const char *skip_blank_db(const char *p) {
      while (*p == ' ' || *p == '\t') ++p;
      return p;
   }
   inline const char *skip_word_db(const char *p) {
      while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') ++p;
      return p;
   }
   inline bool is_skip_line_db(const char *line) {
      const char *p = skip_blank_db(line);
      return *p == '\0' || *p == '#' || *p == '\n' || *p == '\r';
   }

   // Read every numbered N.table file in a directory; returns sorted ints.
   // (Replaces the gemmi original's popen("ls -1 ...") file enumeration.)
   std::vector<int> list_table_files(const std::string &dir) {
      std::vector<int> out;
      if (!std::filesystem::is_directory(dir))
         return out;
      for (const auto &entry : std::filesystem::directory_iterator(dir)) {
         if (!entry.is_regular_file()) continue;
         int n = 0;
         std::string name = entry.path().filename().string();
         if (std::sscanf(name.c_str(), "%d.table", &n) == 1)
            out.push_back(n);
      }
      std::sort(out.begin(), out.end());
      return out;
   }

   // Load the coded -> full-type atom-type map (allAtomTypesFromMolsCoded.list).
   std::unordered_map<std::string, std::string>
   load_atom_codes(const std::string &path) {
      std::unordered_map<std::string, std::string> out;
      FILE *f = std::fopen(path.c_str(), "r");
      if (!f) return out;
      char line[2048];
      while (std::fgets(line, sizeof(line), f)) {
         if (is_skip_line_db(line)) continue;
         const char *p = line;
         const char *s = skip_blank_db(p); p = skip_word_db(s);
         std::string code(s, p - s);
         s = skip_blank_db(p); p = skip_word_db(s);
         std::string full(s, p - s);
         if (!code.empty() && !full.empty())
            out.emplace(std::move(code), std::move(full));
      }
      std::fclose(f);
      return out;
   }

   inline std::string prefix_before_db(const std::string &s, char c) {
      auto pos = s.find(c);
      return pos == std::string::npos ? s : s.substr(0, pos);
   }

   // --- bond table converter -----------------------------------------------

   void convert_bond_tables(SqliteDB &db, const std::string &bond_dir,
                            const std::unordered_map<std::string, std::string> &codes) {
      db.exec(
         "DROP TABLE IF EXISTS bond_entries;"
         "CREATE TABLE bond_entries ("
         "  ha1 INTEGER, ha2 INTEGER,"
         "  hybr_comb TEXT, in_ring TEXT,"
         "  a1_nb2 TEXT, a2_nb2 TEXT, a1_nb TEXT, a2_nb TEXT,"
         "  a1_type_m TEXT, a2_type_m TEXT,"
         "  a1_type_f TEXT, a2_type_f TEXT,"
         "  value REAL, sigma REAL, count INTEGER,"
         "  value_1d REAL, sigma_1d REAL, count_1d INTEGER"
         ");"
      );

      SqliteStmt ins;
      ins.prepare(db.db,
         "INSERT INTO bond_entries (ha1, ha2, hybr_comb, in_ring,"
         "  a1_nb2, a2_nb2, a1_nb, a2_nb,"
         "  a1_type_m, a2_type_m, a1_type_f, a2_type_f,"
         "  value, sigma, count, value_1d, sigma_1d, count_1d)"
         " VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
      );

      db.exec("BEGIN");

      int n_files = 0, n_rows = 0;
      for (int file_num : list_table_files(bond_dir)) {
         std::string path = bond_dir + "/" + std::to_string(file_num) + ".table";
         FILE *f = std::fopen(path.c_str(), "r");
         if (!f) continue;
         ++n_files;

         char line[512];
         while (std::fgets(line, sizeof(line), f)) {
            if (is_skip_line_db(line)) continue;

            const char *p = line;
            const char *s;

            // ha1, ha2
            int ha1 = std::atoi(p);
            while (*p && *p != ' ' && *p != '\t') ++p;
            p = skip_blank_db(p);
            int ha2 = std::atoi(p);
            while (*p && *p != ' ' && *p != '\t') ++p;

            s = skip_blank_db(p); p = skip_word_db(s); std::string hybr_comb(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string in_ring(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a1_nb2(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a2_nb2(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a1_nb(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a2_nb(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string code1(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string code2(s, p - s);
            if (code2.empty()) continue;

            char *endp = nullptr;
            double value   = std::strtod(p, &endp); p = endp;
            double sigma   = std::strtod(p, &endp); p = endp;
            int    count   = std::strtol(p, &endp, 10); p = endp;
            double value2  = std::strtod(p, &endp); p = endp;
            double sigma2  = std::strtod(p, &endp); p = endp;
            int    count2  = std::strtol(p, &endp, 10); p = endp;

            auto it1 = codes.find(code1);
            auto it2 = codes.find(code2);
            std::string a1_type_f = it1 != codes.end() ? it1->second : std::string();
            std::string a2_type_f = it2 != codes.end() ? it2->second : std::string();
            std::string a1_type_m = prefix_before_db(a1_type_f, '{');
            std::string a2_type_m = prefix_before_db(a2_type_f, '{');

            sqlite3_reset(ins.st);
            sqlite3_bind_int   (ins.st, 1,  ha1);
            sqlite3_bind_int   (ins.st, 2,  ha2);
            sqlite3_bind_text  (ins.st, 3,  hybr_comb.c_str(), -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, 4,  in_ring.c_str(),   -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, 5,  a1_nb2.c_str(),    -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, 6,  a2_nb2.c_str(),    -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, 7,  a1_nb.c_str(),     -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, 8,  a2_nb.c_str(),     -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, 9,  a1_type_m.c_str(), -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, 10, a2_type_m.c_str(), -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, 11, a1_type_f.c_str(), -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, 12, a2_type_f.c_str(), -1, SQLITE_TRANSIENT);
            sqlite3_bind_double(ins.st, 13, value);
            sqlite3_bind_double(ins.st, 14, sigma);
            sqlite3_bind_int   (ins.st, 15, count);
            sqlite3_bind_double(ins.st, 16, value2);
            sqlite3_bind_double(ins.st, 17, sigma2);
            sqlite3_bind_int   (ins.st, 18, count2);
            if (sqlite3_step(ins.st) != SQLITE_DONE)
               throw std::runtime_error(std::string("acedrg-db: bond insert failed: ")
                                        + sqlite3_errmsg(db.db));
            ++n_rows;
         }
         std::fclose(f);
      }

      db.exec("COMMIT");
      std::cout << "    bond_entries: " << n_files << " files, " << n_rows << " rows" << std::endl;

      // Build indices last (faster than maintaining them during bulk insert).
      db.exec("CREATE INDEX idx_bond_hash ON bond_entries (ha1, ha2);");
      db.exec("CREATE INDEX idx_bond_full ON bond_entries"
              " (ha1, ha2, hybr_comb, in_ring,"
              "  a1_nb2, a2_nb2, a1_nb, a2_nb);");
   }

   // --- angle table converter -----------------------------------------------

   void convert_angle_tables(SqliteDB &db, const std::string &angle_dir,
                             const std::unordered_map<std::string, std::string> &codes) {
      db.exec(
         "DROP TABLE IF EXISTS angle_entries;"
         "CREATE TABLE angle_entries ("
         "  ha1 INTEGER, ha2 INTEGER, ha3 INTEGER,"
         "  value_key TEXT,"
         "  a1_root TEXT, a2_root TEXT, a3_root TEXT,"
         "  a1_nb2 TEXT, a2_nb2 TEXT, a3_nb2 TEXT,"
         "  a1_nb  TEXT, a2_nb  TEXT, a3_nb  TEXT,"
         "  a1_type TEXT, a2_type TEXT, a3_type TEXT,"
         "  v1 REAL, s1 REAL, c1 INTEGER,"
         "  v2 REAL, s2 REAL, c2 INTEGER,"
         "  v3 REAL, s3 REAL, c3 INTEGER,"
         "  v4 REAL, s4 REAL, c4 INTEGER,"
         "  v5 REAL, s5 REAL, c5 INTEGER,"
         "  v6 REAL, s6 REAL, c6 INTEGER"
         ");"
      );

      SqliteStmt ins;
      ins.prepare(db.db,
         "INSERT INTO angle_entries VALUES ("
         "?,?,?,?,"             // ha1..3, value_key
         "?,?,?,?,?,?,?,?,?,"   // 3x (root, nb2, nb)
         "?,?,?,"               // types
         "?,?,?,?,?,?,?,?,?,"   // 3 levels of v/s/c
         "?,?,?,?,?,?,?,?,?)"   // 3 more levels
      );

      db.exec("BEGIN");

      int n_files = 0, n_rows = 0;
      for (int file_num : list_table_files(angle_dir)) {
         std::string path = angle_dir + "/" + std::to_string(file_num) + ".table";
         FILE *f = std::fopen(path.c_str(), "r");
         if (!f) continue;
         ++n_files;

         char line[1024];
         while (std::fgets(line, sizeof(line), f)) {
            if (is_skip_line_db(line)) continue;

            const char *p = line;
            const char *s;
            int ha1 = std::atoi(p); while (*p && *p != ' ' && *p != '\t') ++p; p = skip_blank_db(p);
            int ha2 = std::atoi(p); while (*p && *p != ' ' && *p != '\t') ++p; p = skip_blank_db(p);
            int ha3 = std::atoi(p); while (*p && *p != ' ' && *p != '\t') ++p;

            s = skip_blank_db(p); p = skip_word_db(s); std::string value_key(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a1_root(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a2_root(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a3_root(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a1_nb2(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a2_nb2(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a3_nb2(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a1_nb(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a2_nb(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string a3_nb(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string code1(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string code2(s, p - s);
            s = skip_blank_db(p); p = skip_word_db(s); std::string code3(s, p - s);
            if (code3.empty()) continue;

            double v[6]; double sg[6]; int c[6];
            char *endp = nullptr;
            for (int i = 0; i < 6; ++i) {
               v[i]  = std::strtod(p, &endp); p = endp;
               sg[i] = std::strtod(p, &endp); p = endp;
               c[i]  = std::strtol(p, &endp, 10); p = endp;
            }

            auto it1 = codes.find(code1);
            auto it2 = codes.find(code2);
            auto it3 = codes.find(code3);
            std::string a1_type = it1 != codes.end() ? prefix_before_db(it1->second, '{') : std::string();
            std::string a2_type = it2 != codes.end() ? prefix_before_db(it2->second, '{') : std::string();
            std::string a3_type = it3 != codes.end() ? prefix_before_db(it3->second, '{') : std::string();

            sqlite3_reset(ins.st);
            int b = 1;
            sqlite3_bind_int   (ins.st, b++, ha1);
            sqlite3_bind_int   (ins.st, b++, ha2);
            sqlite3_bind_int   (ins.st, b++, ha3);
            sqlite3_bind_text  (ins.st, b++, value_key.c_str(), -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a1_root.c_str(),   -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a2_root.c_str(),   -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a3_root.c_str(),   -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a1_nb2.c_str(),    -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a2_nb2.c_str(),    -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a3_nb2.c_str(),    -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a1_nb.c_str(),     -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a2_nb.c_str(),     -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a3_nb.c_str(),     -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a1_type.c_str(),   -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a2_type.c_str(),   -1, SQLITE_TRANSIENT);
            sqlite3_bind_text  (ins.st, b++, a3_type.c_str(),   -1, SQLITE_TRANSIENT);
            for (int i = 0; i < 6; ++i) {
               sqlite3_bind_double(ins.st, b++, v[i]);
               sqlite3_bind_double(ins.st, b++, sg[i]);
               sqlite3_bind_int   (ins.st, b++, c[i]);
            }
            if (sqlite3_step(ins.st) != SQLITE_DONE)
               throw std::runtime_error(std::string("acedrg-db: angle insert failed: ")
                                        + sqlite3_errmsg(db.db));
            ++n_rows;
         }
         std::fclose(f);
      }

      db.exec("COMMIT");
      std::cout << "    angle_entries: " << n_files << " files, " << n_rows << " rows" << std::endl;

      db.exec("CREATE INDEX idx_angle_hash ON angle_entries (ha1, ha2, ha3, value_key);");
   }

   // Copy the cheap ASCII table files -- everything at the top level of
   // acedrg_ascii_tables_dir except *.sqlite, *.bin and
   // allAtomTypesFromMolsCoded.list -- into output_data_dir. Sub-directories
   // (allOrgBondTables/, allOrgAngleTables/, ...) are skipped: they are
   // replaced by acedrg.sqlite and nothing else reads them.
   void copy_cheap_tables(const std::string &acedrg_ascii_tables_dir,
                          const std::string &output_data_dir) {

      int n_copied = 0;
      for (const auto &entry : std::filesystem::directory_iterator(acedrg_ascii_tables_dir)) {
         if (!entry.is_regular_file()) continue;
         std::string name = entry.path().filename().string();
         if (name == "allAtomTypesFromMolsCoded.list") continue;
         if (entry.path().extension() == ".sqlite") continue;
         if (entry.path().extension() == ".bin") continue;
         std::filesystem::path dest = std::filesystem::path(output_data_dir) / name;
         std::filesystem::copy_file(entry.path(), dest,
                                    std::filesystem::copy_options::overwrite_existing);
         ++n_copied;
      }
      std::cout << "    copied " << n_copied << " cheap table files" << std::endl;
   }

} // namespace

bool
coot::acedrg_sqlite_tables::build(const std::string &acedrg_ascii_tables_dir,
                                  const std::string &output_data_dir) {

   try {
      std::filesystem::create_directories(output_data_dir);

      copy_cheap_tables(acedrg_ascii_tables_dir, output_data_dir);

      std::string sqlite_path = output_data_dir + "/acedrg.sqlite";

      // Erase any prior file so we don't append to old contents.
      std::remove(sqlite_path.c_str());

      SqliteDB db;
      db.open(sqlite_path, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);

      // Pragmas to speed bulk insert.
      db.exec("PRAGMA journal_mode = OFF;");
      db.exec("PRAGMA synchronous  = OFF;");
      db.exec("PRAGMA temp_store   = MEMORY;");
      db.exec("PRAGMA cache_size   = -200000;");  // 200 MB page cache during build

      auto codes = load_atom_codes(acedrg_ascii_tables_dir + "/allAtomTypesFromMolsCoded.list");
      if (codes.empty()) {
         std::cout << "coot-make-acedrg-sqlite: no atom-type codes found in "
                   << acedrg_ascii_tables_dir << "/allAtomTypesFromMolsCoded.list" << std::endl;
         return false;
      }
      std::cout << "    atom-type codes: " << codes.size() << " entries" << std::endl;

      convert_bond_tables(db, acedrg_ascii_tables_dir + "/allOrgBondTables", codes);
      convert_angle_tables(db, acedrg_ascii_tables_dir + "/allOrgAngleTables", codes);

      // Final analysis pass so the query planner has stats from the start.
      db.exec("ANALYZE;");

      return true;
   }
   catch (const std::exception &e) {
      std::cout << "coot-make-acedrg-sqlite: " << e.what() << std::endl;
      return false;
   }
}

#else // !USE_SQLITE3

bool
coot::acedrg_sqlite_tables::build(const std::string &acedrg_ascii_tables_dir,
                                  const std::string &output_data_dir) {

   std::cout << "coot-make-acedrg-sqlite: compiled without SQLite3 support" << std::endl;
   return false;
}

#endif // USE_SQLITE3
