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

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif


#include <map>

#include "mmdb_manager.h"
#include "clipper/core/coords.h"
#include "coot-coord-utils.hh"
#include "protein-geometry.hh"
#include "lbg-shared.hh"

#include "flev-annotations.hh"

namespace coot { 

   void write_solvent_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
				      CResidue *reference_residue);

   
   std::map<std::string, std::string> make_flat_ligand_name_map(CResidue *flat_res);

   // return 100 if no other contact found (strange!)
   // 
   double find_water_protein_length(CResidue *ligand_residue, CMMDBManager *mol);


   std::vector<fle_ligand_bond_t> get_covalent_bonds(CMMDBManager *mol,
						     int SelHnd_lig,
						     int SelHnd_all,
						     const residue_spec_t &ligand_spec,
						     const protein_geometry &geom);

   std::vector<fle_ligand_bond_t> get_metal_bonds(CResidue *ligand_res,
						  const std::vector<CResidue *> &residues);


   std::vector<fle_ligand_bond_t> get_fle_ligand_bonds(CResidue *res_ref,
						       const std::vector<CResidue *> &residues,
						       const std::map<std::string, std::string> &name_map);
   // uses the coot::h_bond class (which uses the dictionary).
   // 
   std::vector<fle_ligand_bond_t> get_fle_ligand_bonds(CResidue *res_ref,
						       const std::vector<CResidue *> &residues,
						       CMMDBManager *mol,
						       const std::map<std::string, std::string> &name_map,
						       const protein_geometry &geom);
   
   void write_fle_centres(const std::vector<fle_residues_helper_t> &v,
			  const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
			  const std::vector<coot::solvent_exposure_difference_helper_t> &sed,
			  const pi_stacking_container_t &stacking,
			  CResidue *res_flat);

   bool standard_residue_name_p(const std::string &rn);

   enum { SP_HYBRIDIZATION, SP2_HYBRIDIZATION, SP3_HYBRIDIZATION };

   std::vector<std::pair<std::string, int> >
   get_prodrg_hybridizations(const std::string &prodrg_log_file_name);

   std::vector<std::pair<CAtom *, std::vector<clipper::Coord_orth> > >
   get_cannonball_vectors(CResidue *ligand_res_3d,
			  const coot::dictionary_residue_restraints_t &monomer_restraints);

   enum { H_IS_RIDING, H_IS_ROTATABLE }; // shared between named_torsion_t and flev_attached_hydrogens_t.


   // the reference_residue is the flat residue
   void write_ligand_atom_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
					  const coot::flev_attached_hydrogens_t &attached_hydrogens,
					  CResidue *reference_residue);

} // namespace coot

#endif // C_INTERFACE_LIGANDS_HH
