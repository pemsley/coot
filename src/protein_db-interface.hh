/*
 * src/protein_db-interface.hh
 *
 * Copyright 2011 by University of York
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

#include "protein_db/protein_db.h"  // Kevin's

mmdb::Manager *make_mol(const std::vector<ProteinDB::Chain> &chains, const std::string &chain_id, 
		       int first_resno, bool preserve_residue_names);
mmdb::Manager *make_mol(const ProteinDB::Chain &chain, const std::string &chain_id, int first_resno,
		       bool preserve_residue_names);
std::vector<mmdb::Residue *> 
add_chain_to_molecule(const ProteinDB::Chain &chain, const std::string &chain_id, 
		      int first_res_no, bool preserve_residue_names, mmdb::Manager *mol);

void add_cbs_and_os(std::vector<mmdb::Residue *> needs_cb_and_o, 
		    mmdb::Manager *mol);
