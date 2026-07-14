/* lidia-core/cod-types-similarity.cc
 *
 * Copyright 2026 by Medical Research Council
 * Author: Paul Emsley
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

#include <fstream>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <stdexcept>

#include "cod-types-similarity.hh"

// This file holds only the pure (RDKit-free) parts: normalization, the binary
// database loader, and the Jaccard search.  normalized_types_for_rdkit_mol()
// lives in cod-atom-types.cc, next to the atom typer it depends on, so that
// this translation unit can be compiled and tested without RDKit.

// ---------------------------------------------------------------------------
// Atom type normalization (a faithful port of the Python compute-cod-types.py
// _extract_elements / _parse_neighbor_groups / _normalize_atom_type at
// --level normalized).  Keep these in lock-step with that script.
// ---------------------------------------------------------------------------

namespace {

   // Extract element symbols from a string like "C[6a]C[6a]H" -> {"C","C","H"},
   // skipping over [...] ring annotations.
   std::vector<std::string> extract_elements(const std::string &s) {
      std::vector<std::string> elements;
      std::size_t i = 0;
      while (i < s.size()) {
         char c = s[i];
         if (c == '[') {
            std::size_t j = s.find(']', i);
            if (j == std::string::npos) break;
            i = j + 1;
         } else if (std::isupper(static_cast<unsigned char>(c))) {
            std::string elem(1, c);
            i++;
            if (i < s.size() && std::islower(static_cast<unsigned char>(s[i])) && s[i] != '[') {
               elem += s[i];
               i++;
            }
            elements.push_back(elem);
         } else {
            i++;
         }
      }
      return elements;
   }

   // Match the central atom (element + optional [ring]) at the start of an
   // environment-stripped type string.  Returns the central atom string and
   // sets pos to the first character after it.  Returns "" if no match.
   std::string parse_central(const std::string &atype, std::size_t &pos) {
      pos = 0;
      if (atype.empty() || !std::isupper(static_cast<unsigned char>(atype[0])))
         return std::string();
      std::size_t i = 0;
      std::string central(1, atype[i]);
      i++;
      if (i < atype.size() && std::islower(static_cast<unsigned char>(atype[i]))) {
         central += atype[i];
         i++;
      }
      if (i < atype.size() && atype[i] == '[') {
         std::size_t j = atype.find(']', i);
         if (j != std::string::npos) {
            central += atype.substr(i, j - i + 1);
            i = j + 1;
         }
      }
      pos = i;
      return central;
   }
}

std::string
cod::normalize_atom_type(const std::string &cod_type) {

   // Strip the {...} environment.
   std::string atype = cod_type;
   std::size_t brace = atype.find('{');
   if (brace != std::string::npos)
      atype = atype.substr(0, brace);

   std::size_t pos = 0;
   std::string central = parse_central(atype, pos);
   if (central.empty())
      return atype; // no central match - return the env-stripped string as-is

   const std::string rest = atype.substr(pos);

   // Parse parenthesized neighbour groups with optional trailing repeat count.
   std::vector<std::pair<std::string, int> > groups;
   std::size_t i = 0;
   while (i < rest.size()) {
      if (rest[i] == '(') {
         int depth = 1;
         std::size_t j = i + 1;
         while (j < rest.size() && depth > 0) {
            if (rest[j] == '(') depth++;
            else if (rest[j] == ')') depth--;
            j++;
         }
         // group content is rest[i+1 .. j-2]
         std::string group_content = rest.substr(i + 1, (j - 1) - (i + 1));
         int count = 1;
         std::size_t k = j;
         while (k < rest.size() && std::isdigit(static_cast<unsigned char>(rest[k]))) k++;
         if (k > j)
            count = std::stoi(rest.substr(j, k - j));
         std::vector<std::string> elems = extract_elements(group_content);
         std::sort(elems.begin(), elems.end());
         std::string sorted_elems;
         for (const auto &e : elems) sorted_elems += e;
         groups.push_back(std::make_pair(sorted_elems, count));
         i = k;
      } else {
         i++;
      }
   }

   if (groups.empty())
      return central;

   // Canonical order: sort by (elements, count). std::pair<string,int> does this.
   std::sort(groups.begin(), groups.end());

   std::string out = central;
   for (const auto &g : groups) {
      out += "(";
      out += g.first;
      out += ")";
      if (g.second > 1)
         out += std::to_string(g.second);
   }
   return out;
}

// ---------------------------------------------------------------------------
// Binary database (ACDR format written by compute-cod-types.py:save_db_binary)
//   magic(4) "ACDR"
//   uint32 n_types, uint32 n_molecules                 (little-endian)
//   dictionary:  per type   uint16 len + utf8 bytes
//   molecules:   per mol    uint8 cid_len + cid + uint16 n_ids + n_ids*uint16
// ---------------------------------------------------------------------------

cod::cod_types_db_t
cod::load_cod_types_db_binary(const std::string &path) {

   cod_types_db_t db;
   std::ifstream f(path.c_str(), std::ios::binary);
   if (!f)
      throw std::runtime_error("cannot open COD types DB: " + path);

   char magic[4] = {0,0,0,0};
   f.read(magic, 4);
   if (!f || std::string(magic, 4) != "ACDR")
      throw std::runtime_error("bad magic in COD types DB: " + path);

   auto read_u32 = [&f]() -> std::uint32_t {
      unsigned char b[4] = {0,0,0,0};
      f.read(reinterpret_cast<char *>(b), 4);
      return  static_cast<std::uint32_t>(b[0])        |
             (static_cast<std::uint32_t>(b[1]) <<  8) |
             (static_cast<std::uint32_t>(b[2]) << 16) |
             (static_cast<std::uint32_t>(b[3]) << 24);
   };
   auto read_u16 = [&f]() -> std::uint16_t {
      unsigned char b[2] = {0,0};
      f.read(reinterpret_cast<char *>(b), 2);
      return static_cast<std::uint16_t>(b[0] | (b[1] << 8));
   };
   auto read_u8 = [&f]() -> std::uint8_t {
      unsigned char b = 0;
      f.read(reinterpret_cast<char *>(&b), 1);
      return b;
   };

   std::uint32_t n_types = read_u32();
   std::uint32_t n_mol   = read_u32();

   std::vector<std::string> dict;
   dict.reserve(n_types);
   for (std::uint32_t i = 0; i < n_types; i++) {
      std::uint16_t len = read_u16();
      std::string s(len, '\0');
      if (len) f.read(&s[0], len);
      dict.push_back(s);
   }
   if (!f)
      throw std::runtime_error("truncated COD types DB (dictionary): " + path);

   for (std::uint32_t i = 0; i < n_mol; i++) {
      std::uint8_t cid_len = read_u8();
      std::string cid(cid_len, '\0');
      if (cid_len) f.read(&cid[0], cid_len);
      std::uint16_t n_ids = read_u16();
      std::set<std::string> types;
      for (std::uint16_t k = 0; k < n_ids; k++) {
         std::uint16_t tid = read_u16();
         if (tid < dict.size())
            types.insert(dict[tid]);
      }
      if (!f)
         throw std::runtime_error("truncated COD types DB (molecules): " + path);
      db.molecules[cid] = types;
   }

   return db;
}

// ---------------------------------------------------------------------------
// Jaccard search
// ---------------------------------------------------------------------------

std::vector<cod::cod_match_t>
cod::search(const std::set<std::string> &query_types,
            const cod_types_db_t &db,
            double min_score,
            double beta) {

   std::vector<cod_match_t> results;
   if (query_types.empty())
      return results;

   const int q_size = static_cast<int>(query_types.size());

   for (const auto &m : db.molecules) {
      const std::set<std::string> &ref = m.second;

      // intersection: iterate the smaller set, probe the larger
      const std::set<std::string> &small = (query_types.size() < ref.size()) ? query_types : ref;
      const std::set<std::string> &large = (query_types.size() < ref.size()) ? ref : query_types;
      int common = 0;
      for (const auto &s : small)
         if (large.count(s)) common++;
      if (common == 0)
         continue;

      // Asymmetric Tversky index, alpha = 1:
      //   score = common / ( common + (|A|-common) + beta*(|B|-common) )
      //         = common / ( |A| + beta*(|B|-common) )
      const int r_size = static_cast<int>(ref.size());
      double denom = static_cast<double>(q_size) + beta * static_cast<double>(r_size - common);
      double s = (denom > 0.0) ? (static_cast<double>(common) / denom) : 0.0;
      if (s > min_score)
         results.push_back(cod_match_t(m.first, s, common));
   }

   std::sort(results.begin(), results.end(),
             [](const cod_match_t &a, const cod_match_t &b) {
                if (a.score != b.score) return a.score > b.score;
                return a.comp_id < b.comp_id;
             });
   return results;
}
