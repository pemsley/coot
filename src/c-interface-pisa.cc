/* src/main.cc
 * 
 * Copyright 2009, 2010 The University of Oxford
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


#include "graphics-info.h"
#include <gtk/gtk.h>
#include "c-interface.h"
#include "mmdb_manager.h"

/*  ----------------------------------------------------------------------- */
/*                  PISA Interface                                      */
/*  ----------------------------------------------------------------------- */

// return the new model_number or -1; 
int pisa_interaction(int imol_1, int imol_2) {

   int imodel_new = -1;

   float dist = 4.0;
   if (is_valid_model_molecule(imol_1)) { 
      if (is_valid_model_molecule(imol_2)) {

	 CMMDBManager *mol1 = graphics_info_t::molecules[imol_1].atom_sel.mol;
	 CMMDBManager *mol2 = graphics_info_t::molecules[imol_2].atom_sel.mol;

	 coot::close_residues_from_different_molecules_t cr;
	 std::pair<std::vector<CResidue *>, std::vector<CResidue *> > res_pair = 
	    cr.close_residues(mol1, mol2, dist);

	 if (res_pair.first.size() > 0) { 
	    std::pair<bool, CMMDBManager *> nm =
	       coot::util::create_mmdbmanager_from_residue_vector(res_pair.first);
	    if (nm.second) {
	       int imol = graphics_info_t::create_molecule();
	       atom_selection_container_t asc = make_asc(nm.second);
	       std::string name = "interacting residues from ";
	       name += coot::util::int_to_string(imol_1);
	       graphics_info_t::molecules[imol].install_model(imol, asc, name, 1);
	       imodel_new = imol;
	    } else {
	       std::cout << "WARNING:: no molecule from create_mmdbmanager_from_residue_vector"
			 << std::endl;
	    } 
	 }
	 
	 if (res_pair.second.size() > 0) { 
	    std::pair<bool, CMMDBManager *> nm =
	       coot::util::create_mmdbmanager_from_residue_vector(res_pair.second);
	    if (nm.second) {
	       int imol = graphics_info_t::create_molecule();
	       atom_selection_container_t asc = make_asc(nm.second);
	       std::string name = "interacting residues from ";
	       name += coot::util::int_to_string(imol_2);
	       graphics_info_t::molecules[imol].install_model(imol, asc, name, 1);
	    } else {
	       std::cout << "WARNING:: no molecule from create_mmdbmanager_from_residue_vector"
			 << std::endl;
	    } 
	 }

	 cr.clean_up();
	 graphics_draw();
      }
   }
   return imodel_new;
} 



int pisa_interface_scm(int imol_1, int imol_2, SCM interface_description_scm) {

   // coot::pisa_molecule_t pisa_molecule_1;
   // coot::pisa_molecule_t pisa_molecule_2;

   // coot::pisa_interface_t pi(imol_1, imol_2, pisa_molecule_1, pisa_molecule_2);
   return -1;
}


void pisa_clear_interfaces() {
   

}


