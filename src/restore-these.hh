/*
 * src/restore-these.hh
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

bool residue_to_sdf_file(int imol, const char *chain_id, int res_no, const char *ins_code, 
			 const char *sdf_file_name, bool kekulize);

bool residue_to_mdl_file_for_mogul(int imol, const char *chain_id,
				   int res_no, const char *ins_code, 
				   const char *mdl_file_name);

bool show_feats(int imol, const char *chain_id, int resno, const char *ins_code);

// rdkit chemical features.
// now wrappes generate_meshes()
bool show_feats(int imol, mmdb::Residue *residue_p, const coot::protein_geometry &geom);

