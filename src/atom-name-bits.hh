/* src/atom-name-bits.cc
 *
 * Copyright 2008 by the University of Oxford
 * Copyright 2015 by Medical Research Council
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef ATOM_NAME_BITS_HH
#define ATOM_NAME_BITS_HH

#include <string>
#include <mmdb2/mmdb_manager.h>
#include "utils/coot-utils.hh"

namespace coot {
   // e.g.
   // "Mg" -> <"  MG", "MG">
   // "I"  -> <"   I", "IOD">
   //
   class atom_name_bits_t {
   public:
      atom_name_bits_t() { filled = false; }
      bool filled;
      std::string atom_name;
      std::string element_name;
      std::string res_name;
      atom_name_bits_t(const std::string &type) {
         filled = false;
         if (type == "Br") {
            atom_name = "BR  ";
            element_name = "BR";
            res_name = "BR";
            filled = true;
         }
         if (type == "Ca") {
            atom_name = "CA  ";
            element_name = "CA";
            res_name = "CA";
            filled = true;
         }
         if (type == "Na") {
            atom_name = "NA  ";
            element_name = "NA";
            res_name = "NA";
            filled = true;
         }
         if (type == "K") {
            atom_name = " K  ";
            element_name = "K";
            res_name = "K";
            filled = true;
         }
         if (type == "Cl") {
            atom_name = "CL  ";
            element_name = "CL";
            res_name = "CL";
            filled = true;
         }
         if (type == "I") {
            atom_name = " I  ";
            element_name = "I";
            res_name = "IOD";
            filled = true;
         }
         if (type == "Mg") {
            atom_name = "MG  ";
            element_name = "MG";
            res_name = "MG";
            filled = true;
         }
         if (type == "Zn" || type == "ZN") {
            atom_name = "ZN  ";
            element_name = "ZN";
            res_name = "ZN";
            filled = true;
         }
         if (! filled) {
            // make up (guess) the residue type and element
            std::string at_name = util::upcase(type);
            atom_name = at_name;
            res_name = at_name;
            if (atom_name.length() == 2)
               atom_name = atom_name + "  ";
            // if (atom_name.length() == 1)
            // atom_name = atom_name + "   ";
            if (atom_name.length() == 1)
               atom_name = std::string(" ") + atom_name + std::string("  ");
            element_name = at_name;
            if (type.length() > 4)
               atom_name = at_name.substr(0,4);
            if (type.length() > 3)
               res_name = at_name.substr(0,3);
            if (type.length() > 2)
               element_name = at_name.substr(0,2);
            filled = true;
         }
      }
      void SetAtom(mmdb::Atom *at, mmdb::Residue *res) {
         if (filled) {
            at->SetAtomName(atom_name.c_str());
            at->SetElementName(element_name.c_str());
            if (false) // can't fix the MMDB bug, I think.
               if (res_name.length() == 2)
                  res_name = res_name + " ";
            res->SetResName(res_name.c_str());
         }
      }
   };
}

#endif // ATOM_NAME_BITS_HH

