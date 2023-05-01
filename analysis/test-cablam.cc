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
