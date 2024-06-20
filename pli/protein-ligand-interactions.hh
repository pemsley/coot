/*
 * pli/protein-ligand-interactions.hh
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

// #include <mmdb2/mmdb_manager.h>

#include "flev-annotations.hh"

namespace pli {

   std::vector<fle_ligand_bond_t>
   protein_ligand_interactions(mmdb::Residue *residue_p,
			       mmdb::Manager *mol,
			       coot::protein_geometry *geom_p,
			       float h_bond_dist_max);


   std::vector<fle_ligand_bond_t> get_covalent_bonds(mmdb::Manager *mol,
						     int SelHnd_lig,
						     int SelHnd_all,
						     const coot::residue_spec_t &ligand_spec,
						     const coot::protein_geometry &geom);
   // which calls 
   std::vector<fle_ligand_bond_t> get_covalent_bonds_by_distance(mmdb::Manager *mol,
						     int SelHnd_lig,
						     int SelHnd_all,
						     const coot::residue_spec_t &ligand_spec,
						     const coot::protein_geometry &geom);
   std::vector<fle_ligand_bond_t> get_covalent_bonds_by_links(mmdb::Residue *residue_ligand_p,
							      mmdb::Manager *mol);

   std::vector<fle_ligand_bond_t> get_metal_bonds(mmdb::Residue *ligand_res,
						  const std::vector<mmdb::Residue *> &residues);


   // uses the coot::h_bond class (which uses the dictionary).
   // 
   std::vector<fle_ligand_bond_t> get_fle_ligand_bonds(mmdb::Residue *res_ref,
						       const std::vector<mmdb::Residue *> &residues,
						       mmdb::Manager *mol,
						       const std::map<std::string, std::string> &name_map,
						       const coot::protein_geometry &geom,
						       float water_dist_max,
						       float h_bond_dist_max);

   // return 100 if no other contact found (strange!)
   // 
   double find_water_protein_length(mmdb::Residue *ligand_residue, mmdb::Manager *mol);

}
