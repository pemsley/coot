/* lidia-core/cod-types-similarity.hh
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

#ifndef COD_TYPES_SIMILARITY_HH
#define COD_TYPES_SIMILARITY_HH

#include <string>
#include <set>
#include <map>
#include <vector>

// Molecule similarity by sets of acedrg/COD atom types.
//
// A molecule is represented by the set of its normalized non-H acedrg/COD
// atom-type strings.  Two molecules are compared by the Jaccard similarity of
// these sets.  The reference database is built offline by compute-cod-types.py;
// the normalization here must reproduce that script's --level normalized output
// exactly, otherwise query types will not match the stored types.

namespace RDKit { class ROMol; } // avoid pulling RDKit headers into this header

namespace cod {

   // Normalize a level-4 acedrg/COD atom-type string to the canonical "normalized"
   // form used by compute-cod-types.py: strip the {...} environment, sort the
   // element letters within each neighbour group, sort the groups, keep the repeat
   // counts, and keep the central atom + [ring] annotation.
   //
   // e.g.  C[6a](C[6a]C[6a]H)2(OC){1|C<3>,2|H<1>}  ->  C[6a](CCH)2(CO)
   std::string normalize_atom_type(const std::string &cod_type);

   // Reference database: comp_id -> set of normalized non-H type strings.
   class cod_types_db_t {
   public:
      std::map<std::string, std::set<std::string> > molecules;
      bool empty() const { return molecules.empty(); }
      std::size_t size() const { return molecules.size(); }
   };

   // Load the ACDR binary produced by compute-cod-types.py.
   // Throws std::runtime_error on a missing or malformed file.
   cod_types_db_t load_cod_types_db_binary(const std::string &path);

   class cod_match_t {
   public:
      std::string comp_id;
      double score;
      int n_common;
      cod_match_t(const std::string &c, double s, int n) : comp_id(c), score(s), n_common(n) {}
   };

   // Compare a query type-set A against every reference B in the database and
   // rank by an asymmetric Tversky index (alpha = 1 fixed):
   //
   //    score = |A n B| / ( |A n B| + |A - B| + beta * |B - A| )
   //
   // alpha = 1 makes missing query types fully penalised, while a small beta
   // (default 0.2) only lightly penalises a reference for being larger than the
   // query - so a known monomer that contains the sketch scores well even when
   // it has many extra atoms, with a mild preference for the tightest fit.
   // (beta = 1 would be Jaccard; beta = 0 would be pure coverage |A n B|/|A|.)
   //
   // Returns the matches with score > min_score, sorted by descending score.
   std::vector<cod_match_t> search(const std::set<std::string> &query_types,
                                   const cod_types_db_t &db,
                                   double min_score,
                                   double beta = 0.2);

   // Compute the normalized non-H type set for an RDKit molecule.  A sanitized
   // copy is typed with cod::atom_types_t, so the input is not modified.
   // Can throw std::runtime_error (from the atom typer).
   std::set<std::string> normalized_types_for_rdkit_mol(const RDKit::ROMol &mol);

}

#endif // COD_TYPES_SIMILARITY_HH
