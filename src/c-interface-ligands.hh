/* src/c-interface-ligands.hh
 * 
 * Copyright 2008, 2009 The University of Oxford
 * Author: Paul Emsley
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
 * write to the Free Software Foundation, Inc., 51 Franklin Street,  02110-1301, USA
 */

#ifndef C_INTERFACE_LIGANDS_HH
#define C_INTERFACE_LIGANDS_HH


#include <map>

#include <mmdb2/mmdb_manager.h>
#include "clipper/core/coords.h"
#include "coot-utils/coot-coord-utils.hh"
#include "geometry/protein-geometry.hh"
#include "lidia-core/lbg-shared.hh"

#include "pli/flev-annotations.hh"
#include "ligand-check.hh"
#include "ideal/simple-restraint.hh"

// not for swig

namespace coot { 

   void write_solvent_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
				      mmdb::Residue *reference_residue);

   
   std::map<std::string, std::string> make_flat_ligand_name_map(mmdb::Residue *flat_res);


   bool standard_residue_name_p(const std::string &rn);

   enum { SP_HYBRIDIZATION, SP2_HYBRIDIZATION, SP3_HYBRIDIZATION };

   std::vector<std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > >
   get_cannonball_vectors(mmdb::Residue *ligand_res_3d,
			  const coot::dictionary_residue_restraints_t &monomer_restraints);


   // the reference_residue is the flat residue
   void write_ligand_atom_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
					  const pli::flev_attached_hydrogens_t &attached_hydrogens,
					  mmdb::Residue *reference_residue);

   void ligand_check_percentiles_dialog(residue_spec_t spec,
					const ligand_report_percentiles_t &lr,
					double percentile_limit);
   // this is a wrapper for the above.
   void ligand_check_dialog(residue_spec_t spec,
			    const ligand_report_absolute_t &lr,
			    double percentile_limit);

   geometry_distortion_info_container_t get_ligand_distortion_summary_info(int imol,
                                                                           residue_spec_t &rs);

} // namespace coot


void coot_contact_dots_for_ligand_internal(int imol, const coot::residue_spec_t &res_spec);

void coot_contact_dots_for_ligand_instancing_version(int imol, const coot::residue_spec_t &res_spec);

// 20230605-PE this was done via scripting previously. Not today though
void smiles_to_simple_3d(const std::string &smiles);

void residue_to_ligand_builder(int imol, const std::string &chain_id, int resno, const std::string &ins_code,
                               double weight_for_3d_distances);
void smiles_to_ligand_builder(const std::string &smiles_string);

#endif // C_INTERFACE_LIGANDS_HH
