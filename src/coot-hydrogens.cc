/*
 * src/coot-hydrogens.cc
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

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/rdkit-interface.hh"
#endif
 
#include <iostream> // fixes undefined strchr, strchrr problems
#include <mmdb2/mmdb_manager.h>

#include "geometry/protein-geometry.hh"
#include "coot-hydrogens.hh"


std::pair<bool, std::string>
coot::add_hydrogens(mmdb::Residue *residue_p,
		    const coot::dictionary_residue_restraints_t &restraints) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   return add_hydrogens_with_rdkit(residue_p, restraints);
#else    
   // return add_hydrogens_with_ccp4_tools(residue_p, restraints);
   return std::pair<bool, std::string> (0, "not implemented");
#endif   

}


