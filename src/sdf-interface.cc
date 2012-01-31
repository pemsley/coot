/* src/lbg-interface.cc
 * 
 * Copyright 2012 by The University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifdef MAKE_ENTERPRISE_TOOLS

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <cstring>
#define ENABLE_NLS // fix dcgettext() header problems on including
		   // libintl.h (via RDKitBase.h etc (including boost
		   // stuff).

#include "graphics-info.h"
#include "sdf-interface.hh"
#include "rdkit-interface.hh"

void residue_to_sdf_file(int imol, const char *chain_id, int res_no, const char *ins_code, 
			 const char *sdf_file_name) {

   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      CResidue *residue_p =
	 graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (residue_p) {
	 try {
	    RDKit::RWMol rdkm = coot::rdkit_mol(residue_p, *g.Geom_p());
	    bool includeStereo = true;
	    int confId = 0;
	    bool kekulize = true;
	    RDKit::MolToMolFile(rdkm, sdf_file_name, includeStereo, confId, kekulize);
	 }
	 catch (std::runtime_error coot_error) {
	    std::cout << coot_error.what() << std::endl;
	    std::string m = "Residue type ";
	    m += residue_p->GetResName();
	    m += " not found in dictionary.";
	 }
	 catch (std::exception rdkit_error) {
	    std::cout << rdkit_error.what() << std::endl;
	 }
      }
   }
} 


#endif 
