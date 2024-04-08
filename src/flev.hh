/*
 * src/flev.hh
 *
 * Copyright 2010 by University of York
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

namespace coot {
   
   std::vector<fle_residues_helper_t>
   get_flev_residue_centres(mmdb::Residue *reference_residue,
			    mmdb::Manager *mol_containing_residue_ligand, 
			    std::vector<mmdb::Residue *> residues,
			    mmdb::Manager *flat_mol);

   // return a vector of the same size as filtered_residues.
   // 
   std::vector<int> make_add_reps_for_near_residues(std::vector<mmdb::Residue *> filtered_residues,
						    int imol);

   void add_animated_ligand_interactions(int imol,
					 const std::vector<fle_ligand_bond_t> &ligand_bonds); 

}


