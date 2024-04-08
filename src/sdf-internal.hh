/*
 * src/sdf-internal.hh
 *
 * Copyright 2012 by University of York
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

#include "lidia-core/rdkit-interface.hh"
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>

namespace chemical_features { 

   // not for public access 
   void show(int imol, const RDKit::ROMol &rdkm, std::string name);
   
   std::pair<bool, clipper::Coord_orth> get_normal_info(RDKit::MolChemicalFeature *feat,
							const RDKit::ROMol &mol,
							const RDKit::Conformer &conf);
   std::pair<bool, clipper::Coord_orth> get_normal_info_aromatic(RDKit::MolChemicalFeature *feat,
								 const RDKit::Conformer &conf);
   std::pair<bool, clipper::Coord_orth> get_normal_info_donor(RDKit::MolChemicalFeature *feat,
							      const RDKit::ROMol &mol,
							      const RDKit::Conformer &conf);
   // return null on inability to make factory.
   RDKit::MolChemicalFeatureFactory *get_feature_factory();
   
}
