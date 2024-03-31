/*
 * analysis/test-cablam.cc
 *
 * Copyright 2013 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <string>
#include "cablam.hh"
#include "typed-distances.hh"

#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"

int main(int argc, char **argv) {

   if (argc > 1) {

      std::string file_name(argv[1]);
      atom_selection_container_t asc = get_atom_selection(file_name, false, false, false);

      if (false) {

	 int n_selected_residues;
	 mmdb::PResidue *SelResidues = 0;
	 int selHnd = asc.mol->NewSelection();
	 asc.mol->Select ( selHnd, mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
			   "*", // Chain id
			   mmdb::ANY_RES,"*",  // starting res
			   mmdb::ANY_RES,"*",  // ending res
			   "*",  // residue name
			   "*",  // Residue must contain this atom name?
			   "*",  // Residue must contain this Element?
			   "*",  // altLocs
			   mmdb::SKEY_NEW // selection key
			   );
	 asc.mol->GetSelIndex (selHnd, SelResidues, n_selected_residues);

	 int imodel = 1;
	 mmdb::Model *model_p = asc.mol->GetModel(imodel);
	 if (model_p)
	    model_p->CalcSecStructure(imodel);

	 coot::cablam c(SelResidues, n_selected_residues);
      }

      if (true) {

	 coot::typed_distances td(asc.mol);
	 td.output();

      }
   }

   return 0;
} 
