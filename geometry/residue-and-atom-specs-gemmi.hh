/* geometry/residue-and-atom-specs-gemmi.hh
 *
 * Copyright 2026 by Paul Emsley
 * Copyright 2026 by Medical Research Council
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

// Gemmi constructors and lookup methods for atom_spec_t and residue_spec_t.
// These extend the existing spec classes with gemmi type support.

#ifndef RESIDUE_AND_ATOM_SPECS_GEMMI_HH
#define RESIDUE_AND_ATOM_SPECS_GEMMI_HH

#include <string>
#include <gemmi/model.hpp>
#include "residue-and-atom-specs.hh"

namespace coot {

   // ==================== atom_spec_t from gemmi ====================

   // Construct an atom_spec_t from gemmi types.
   // Requires the chain and residue context since gemmi::Atom doesn't
   // store its own chain/residue (unlike mmdb::Atom which navigates
   // up the hierarchy).
   atom_spec_t atom_spec_from_gemmi(const gemmi::Chain &chain,
                                    const gemmi::Residue &res,
                                    const gemmi::Atom &at,
                                    int model_number = 1);

   // Does this atom_spec match the given gemmi atom (in context)?
   bool matches_spec(const atom_spec_t &spec,
                     const gemmi::Chain &chain,
                     const gemmi::Residue &res,
                     const gemmi::Atom &at);

   // Find the atom in the structure. Returns a const pointer, or nullptr.
   const gemmi::Atom *get_atom(const atom_spec_t &spec,
                               const gemmi::Structure &st);

   // ==================== residue_spec_t from gemmi ====================

   // Construct a residue_spec_t from gemmi types.
   residue_spec_t residue_spec_from_gemmi(const gemmi::Chain &chain,
                                          const gemmi::Residue &res,
                                          int model_number = 1);

   // Find the residue in the structure. Returns a const pointer, or nullptr.
   const gemmi::Residue *get_residue(const residue_spec_t &spec,
                                     const gemmi::Structure &st);

   // ==================== link_atoms from gemmi ====================

   std::pair<atom_spec_t, atom_spec_t> link_atoms_from_gemmi(const gemmi::Connection &con,
                                                             int model_number = 1);

} // namespace coot

#endif // RESIDUE_AND_ATOM_SPECS_GEMMI_HH
