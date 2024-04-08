/*
 * src/select-atom-info.hh
 *
 * Copyright 2007 by University of York
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

#ifndef SELECT_ATOM_INFO_HH
#define SELECT_ATOM_INFO_HH

#include <mmdb2/mmdb_manager.h>
#include <string>

namespace coot {

   class select_atom_info {
      bool b_factor_editted;
      bool occ_editted;
      bool altloc_edited;
   public:
      int udd;
      int molecule_number;
      std::string chain_id;
      int residue_number;
      std::string insertion_code;
      std::string atom_name;
      std::string altconf;
      float b_factor;
      float occ;
      std::string altloc_new; // the new altloc
      select_atom_info(int udd_in,
                       int molecule_number_in,
                       const std::string &chain_id_in,
                       int residue_number_in,
                       const std::string &insertion_code_in,
                       const std::string &atom_name_in,
                       const std::string &altconf_in) : chain_id(chain_id_in), insertion_code(insertion_code_in),
                                                        atom_name(atom_name_in), altconf(altconf_in) {
        udd = udd_in;
        molecule_number = molecule_number_in;
        residue_number = residue_number_in;
        b_factor = 0.0;
        occ = 0.0;
        //
        occ_editted = false;
        b_factor_editted = false;
        altloc_edited = false;
      }
      select_atom_info() {}

      // add functions that set the values of the edits
      void add_b_factor_edit(float b_factor_in) {
         b_factor_editted = 1;
         b_factor = b_factor_in;
      }
      void add_occ_edit(float occ_in) {
         occ_editted = 1;
         if ((occ_in > 1.0) || (occ_in < 0.0))
            occ = 1.0;
         else
            occ = occ_in;
      }
      void add_altloc_edit(const std::string &s) {
         altloc_edited = 1;
         altloc_new = s;
      }
      bool has_b_factor_edit() const { return b_factor_editted; }
      bool has_occ_edit() const { return occ_editted; }
      bool has_altloc_edit() const { return altloc_edited; }
      // return NULL on atom not found:
      mmdb::Atom *get_atom(mmdb::Manager *mol) const {
         mmdb::Atom *at = NULL;
         if (mol) {
            int SelectionHandle = mol->NewSelection();
            int n_atoms;
            mmdb::PAtom *atoms;
            mol->SelectAtoms(SelectionHandle, 0,
                             chain_id.c_str(),
                             residue_number, insertion_code.c_str(),
                             residue_number, insertion_code.c_str(),
                             "*",
                             atom_name.c_str(), "*",
                             altconf.c_str());
            mol->GetSelIndex(SelectionHandle, atoms, n_atoms);
            if (n_atoms > 0)
               at = atoms[0];
            mol->DeleteSelection(SelectionHandle);
         }
         return at;
      }
   };

}

#endif // SELECT_ATOM_INFO_HH
