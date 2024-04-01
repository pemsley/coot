/*
 * geometry/mol-utils.cc
 *
 * Copyright 2017 by Medical Research Council
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

#include "mol-utils.hh"


std::vector<std::string>
coot::util::get_residue_alt_confs(mmdb::Residue *res) {

   std::vector<std::string> v;
   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   res->GetAtomTable(residue_atoms, nResidueAtoms);
   bool ifound = 0;
   for (int iat=0; iat<nResidueAtoms; iat++) {
      ifound = 0;
      for(unsigned int i=0; i<v.size(); i++) {
         if (std::string(residue_atoms[iat]->altLoc) == v[i]) {
            ifound = 1;
            break;
         }
      }
      if (! ifound)
         v.push_back(std::string(residue_atoms[iat]->altLoc));
   }
   return v;
}


