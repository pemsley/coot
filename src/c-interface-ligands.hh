/* src/c-interface-ligands.hh
 * 
 * Copyright 2008, 2009 The University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#ifndef C_INTERFACE_LIGANDS_HH
#define C_INTERFACE_LIGANDS_HH


#include <map>

#include <mmdb2/mmdb_manager.h>
#include "clipper/core/coords.h"
#include "coot-utils/coot-coord-utils.hh"
#include "geometry/protein-geometry.hh"
#include "lidia-core/lbg-shared.hh"

#include "lbg/flev-annotations.hh"
#include "ligand-check.hh"

namespace coot { 

   void write_solvent_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
				      mmdb::Residue *reference_residue);

   
   std::map<std::string, std::string> make_flat_ligand_name_map(mmdb::Residue *flat_res);

   // return 100 if no other contact found (strange!)
   // 
   double find_water_protein_length(mmdb::Residue *ligand_residue, mmdb::Manager *mol);


   std::vector<fle_ligand_bond_t> get_covalent_bonds(mmdb::Manager *mol,
						     int SelHnd_lig,
						     int SelHnd_all,
						     const residue_spec_t &ligand_spec,
						     const protein_geometry &geom);
   // which calls 
   std::vector<fle_ligand_bond_t> get_covalent_bonds_by_distance(mmdb::Manager *mol,
						     int SelHnd_lig,
						     int SelHnd_all,
						     const residue_spec_t &ligand_spec,
						     const protein_geometry &geom);
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
						       const protein_geometry &geom,
						       float water_dist_max,
						       float h_bond_dist_max);

   bool standard_residue_name_p(const std::string &rn);

   enum { SP_HYBRIDIZATION, SP2_HYBRIDIZATION, SP3_HYBRIDIZATION };

   std::vector<std::pair<std::string, int> >
   get_prodrg_hybridizations(const std::string &prodrg_log_file_name);

   std::vector<std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > >
   get_cannonball_vectors(mmdb::Residue *ligand_res_3d,
			  const coot::dictionary_residue_restraints_t &monomer_restraints);

   enum { H_IS_RIDING, H_IS_ROTATABLE }; // shared between named_torsion_t and flev_attached_hydrogens_t.


   // the reference_residue is the flat residue
   void write_ligand_atom_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
					  const coot::flev_attached_hydrogens_t &attached_hydrogens,
					  mmdb::Residue *reference_residue);

   void ligand_check_percentiles_dialog(residue_spec_t spec,
					const ligand_report_percentiles_t &lr,
					double percentile_limit);
   // this is a wrapper for the above.
   void ligand_check_dialog(residue_spec_t spec,
			    const ligand_report_absolute_t &lr,
			    double percentile_limit);

} // namespace coot

#endif // C_INTERFACE_LIGANDS_HH
