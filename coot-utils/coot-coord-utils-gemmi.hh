/* coot-utils/coot-coord-utils-gemmi.hh
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

// Gemmi equivalents of functions from coot-coord-utils.hh
// These operate on gemmi types instead of mmdb types.

#ifndef COOT_COORD_UTILS_GEMMI_HH
#define COOT_COORD_UTILS_GEMMI_HH

#include <vector>
#include <string>
#include <map>
#include <set>
#include <utility>

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   // Trim padded atom names (and element names) in a gemmi::Structure.
   // Workaround for mmdb padding label_atom_id when reading PDB files.
   // This should eventually be fixed in gemmi's copy_from_mmdb() instead.
   void trim_atom_names(gemmi::Structure &st);

   // ==================== atom-level ====================

   double distance(const gemmi::Atom &at_1, const gemmi::Atom &at_2);
   double angle(const gemmi::Atom &at_1, const gemmi::Atom &at_2, const gemmi::Atom &at_3);
   clipper::Coord_orth co(const gemmi::Atom &at);
   bool is_hydrogen_atom(const gemmi::Atom &at);

   // ==================== residue-level ====================

   namespace util {

      std::pair<bool, clipper::Coord_orth> get_residue_centre(const gemmi::Residue &res);
      std::pair<bool, clipper::Coord_orth> get_CA_position_in_residue(const gemmi::Residue &res);
      std::pair<bool, clipper::Coord_orth> get_CB_position_in_residue(const gemmi::Residue &res);
      short int is_nucleotide(const gemmi::Residue &r);
      bool residue_has_hydrogens_p(const gemmi::Residue &res);
      int residue_has_hetatms(const gemmi::Residue &res);

   } // namespace util

   // ==================== chain-level ====================

   namespace util {

      std::pair<int, int> min_and_max_residues(const gemmi::Chain &chain);
      std::pair<bool, int> min_resno_in_chain(const gemmi::Chain &chain);
      std::pair<bool, int> max_resno_in_chain(const gemmi::Chain &chain);
      std::vector<std::string> residue_types_in_chain(const gemmi::Chain &chain);
      std::pair<unsigned int, unsigned int> get_number_of_protein_or_nucleotides(const gemmi::Chain &chain);

   } // namespace util

   // ==================== structure-level ====================

   std::pair<bool, clipper::Coord_orth> centre_of_molecule(const gemmi::Structure &st);
   std::pair<bool, double> radius_of_gyration(const gemmi::Structure &st);
   bool mol_has_symmetry(const gemmi::Structure &st);
   bool mol_is_anisotropic(const gemmi::Structure &st);
   float get_position_hash(const gemmi::Structure &st);

   namespace util {

      std::vector<std::string> residue_types_in_molecule(const gemmi::Structure &st);
      std::vector<std::string> non_standard_residue_types_in_molecule(const gemmi::Structure &st);
      std::vector<std::string> chains_in_molecule(const gemmi::Structure &st);
      int number_of_residues_in_molecule(const gemmi::Structure &st);
      int max_number_of_residues_in_chain(const gemmi::Structure &st);
      int number_of_chains(const gemmi::Structure &st);
      std::pair<bool, int> max_resno_in_molecule(const gemmi::Structure &st);
      int max_min_max_residue_range(const gemmi::Structure &st);
      std::vector<std::string> alt_confs_in_molecule(const gemmi::Structure &st);
      std::pair<clipper::Coord_orth, clipper::Coord_orth> extents(const gemmi::Structure &st);
      clipper::Coord_orth median_position(const gemmi::Structure &st);
      int count_cis_peptides(const gemmi::Structure &st);

   } // namespace util

} // namespace coot

#endif // COOT_COORD_UTILS_GEMMI_HH
