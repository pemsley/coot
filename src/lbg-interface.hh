/*
 * src/lbg-interface.hh
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

#if 0 // 20230920-PE Now we have layla. Remove this header

// 20140226: now we change things so that the interface functions always get generated,
//           and what the function does depends on MAKE_ENHANCED_LIGAND_TOOLS
// 
// #ifdef MAKE_ENHANCED_LIGAND_TOOLS
void residue_to_ligand_builder(int imol, const char *chain_id, int resno, const char *ins_code, double weight_for_3d_distances);
void smiles_to_ligand_builder(const char *smiles_string);
// #endif

#endif
