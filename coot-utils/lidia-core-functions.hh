/*
 * coot-utils/lidia-core-functions.hh
 *
 * Copyright 2017 by Medical Research Council
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

#ifndef LIDIA_CORE_FUNCTIONS_HH
#define LIDIA_CORE_FUNCTIONS_HH

#include "atom-selection-container.hh"
// #include "lidia-core/lig-build.hh"
#include "lidia-core/lbg-molfile.hh"

namespace coot {
  // mdl mol file support
  atom_selection_container_t mdl_mol_to_asc(const lig_build::molfile_molecule_t &m);
  atom_selection_container_t mdl_mol_to_asc(const lig_build::molfile_molecule_t &m, float b_factor);
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
  atom_selection_container_t mol_to_asc_rdkit(const std::string &file_name);
#endif // MAKE_ENHANCED_LIGAND_TOOLS
}


#endif // LIDIA_CORE_FUNCTIONS_HH

