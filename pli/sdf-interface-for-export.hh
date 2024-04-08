/*
 * pli/sdf-interface-for-export.hh
 *
 * Copyright 2023 by Medical Research Council
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
#ifndef SDF_INTERFACE_FOR_EXPORT_HH
#define SDF_INTERFACE_FOR_EXPORT_HH

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include "lidia-core/rdkit-interface.hh"
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include "coot-utils/simple-mesh.hh"

namespace chemical_features {

   std::vector<coot::simple_mesh_t> generate_meshes(int imol, mmdb::Residue *residue_p, const coot::protein_geometry &geom);

   // choose iconf = 0 if unsure.
   std::vector<coot::simple_mesh_t> generate_meshes(int imol, const RDKit::ROMol &rdkm, int iconf, const std::string &name);
}

#endif // MAKE_ENHANCED_LIGAND_TOOLS

#endif // SDF_INTERFACE_FOR_EXPORT_HH
