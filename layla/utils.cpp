/* layla/utils.cpp
 * 
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
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

#include "utils.hpp"
#include <rdkit/GraphMol/MolOps.h>
#ifndef __EMSCRIPTEN__
#include <glib-2.0/glib.h>
#else
#include "../lhasa/glog_replacement.hpp"
#endif

void coot::layla::remove_non_polar_hydrogens(RDKit::RWMol& mol) {

#if 0 // 20250117-PE
   std::vector<RDKit::Atom*> atoms_to_be_removed;
   auto atoms = mol.atoms();
   for(RDKit::Atom* atom: atoms) {
      if(atom->getAtomicNum() == 1) {
         if(atom->getFormalCharge() == 0) {
            atoms_to_be_removed.push_back(atom);
         }
      }
   }

   for(RDKit::Atom* atom: atoms_to_be_removed) {
      mol.removeAtom(atom);
      try {
         RDKit::MolOps::sanitizeMol(mol);
      } catch (std::exception& e) {
         g_warning("Could not sanitize molecule while removing non-polar hydrogens: %s", e.what());
      }
   }
#endif

   /*
     {
     std::string s = "c1cccnc1";
     RDKit::RWMol *rdkmp = RDKit::SmilesToMol(s);
     RDKit::ROMol rdkm = *rdkmp;
     RDKit::ROMol *rdk_mol_with_no_Hs_ro = RDKit::MolOps::removeHs(rdkm);
     }
   */

   try {
      // so if we pass an ROMol to removeHs(), it returns a value, if we passs a RWMol
      // then it does not.
      //
      bool implicit_only = false;
      bool update_explicit_count = true;
      bool sanitize = true;
      RDKit::MolOps::removeHs(mol, implicit_only, update_explicit_count, sanitize);
   }
   catch (const RDKit::AtomValenceException &e) {
	   std::cout << "WARNING:: " << e.what() << std::endl;
   }

}
