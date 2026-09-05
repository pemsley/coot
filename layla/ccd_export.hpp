/* layla/ccd_export.hpp
 *
 * Copyright 2026 by Global Phasing Ltd.
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

#ifndef LAYLA_CCD_EXPORT_HPP
#define LAYLA_CCD_EXPORT_HPP

#include <string>

namespace RDKit { class ROMol; }

namespace coot::layla {

   // Build a minimal CCD-style mmCIF - the "acedrg -c" input format - from
   // the canvas molecule: kekulized bond orders with pdbx_aromatic_flag,
   // atom names preserved from the molecule's atom "name" properties (a
   // Coot-imported monomer keeps its names; unnamed atoms and the added
   // hydrogens get generated element+counter names), ETKDG-embedded ideal
   // coordinates, and the hand-drawn 2D layout written as the PDBe
   // depiction loops.
   //
   // A C++ port of python/convert-pkl-to-mmcif.py - going via a pickle is
   // not needed when we already hold the RDKit molecule.
   //
   // Returns the mmCIF document as a string. Throws std::runtime_error on
   // failure (no atoms, 3D embedding failure).
   std::string make_acedrg_input_mmcif(const RDKit::ROMol &mol_2d, const std::string &comp_id);

} // namespace coot::layla

#endif // LAYLA_CCD_EXPORT_HPP
