/* src/graphics-info-residues.cc
 * 
 * Copyright 2011 by The University of Oxford.
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
 * Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <algorithm>
#include "graphics-info.h"
#include "interface.h" // for create_multi_residue_torsion_dialog()

void
graphics_info_t:: multi_torsion_residues(int imol, const std::vector<coot::residue_spec_t> &v) {

   if (is_valid_model_molecule(imol)) {

      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (! mol)
	 std::cout << "no (reference) mol" << std::endl;

      // we probably don't need residues, but we do need to check that
      // we have restraints for each of the residue types.
      // 
      std::vector<std::string> residue_types;
      for (unsigned int i=0; i<v.size(); i++) { 
	 CResidue *r= molecules[imol].get_residue(v[i]);
	 if (r) {
	    std::string comp_id = r->GetResName();
	    if (std::find(residue_types.begin(), residue_types.end(), comp_id) ==
		residue_types.end()) {
	       residue_types.push_back(comp_id);
	    }
	 } else {
	    std::cout << "WARNING:: residues " << v[i] << " not found in molecule "
		      << imol << std::endl;
	 }
      }

      // uses dynamic add
      bool dict_status = geom_p->have_dictionary_for_residue_types(residue_types);

      CMMDBManager *moving_mol = coot::util::create_mmdbmanager_from_residue_specs(v, mol);
      if (! moving_mol) {
	 std::cout << "WARNING:: multi_torsion_residues() no moving mol" << std::endl;
      } else {
	 imol_moving_atoms = imol;
	 // now select everything in moving mol
	 int selhnd = moving_mol->NewSelection();
	 PPCAtom atom_selection = 0;
	 int n_selected_atoms;
	 moving_mol->SelectAtoms(selhnd, 0, "*",
				 ANY_RES, "*",
				 ANY_RES, "*",
				 "*", // residue name
				 "*",
				 "*", 
				 "*"); // alt-loc
	 moving_mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);

	 try { 
	    std::vector<std::pair<CAtom *, CAtom *> > pairs = 
	       coot::torsionable_bonds(mol, atom_selection, n_selected_atoms, Geom_p());
	    
	    GtkWidget *w = wrapped_create_multi_residue_torsion_dialog(pairs);
	    gtk_widget_show(w);

	    moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
	    make_moving_atoms_graphics_object(make_asc(moving_mol));
	 }
	 catch (std::runtime_error rte) {
	    std::cout << "WARNING:: " << rte.what() << std::endl;
	 } 

	 moving_mol->DeleteSelection(selhnd);
      }
   }
} 

