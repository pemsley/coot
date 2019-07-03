/* coords/mmdb.h
 * 
 * Copyright 2005 by The University of York
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

#ifndef MMDB_H
#define MMDB_H

#include "lidia-core/lbg-molfile.hh"

// needs to have included mmdb_manager.h"
// and "mmdb-extras.h" for atom_selection_container_t
#include "Cartesian.h"


coot::Cartesian
centre_of_molecule(atom_selection_container_t SelAtom); 



int write_atom_selection_file(atom_selection_container_t asc,
			      const std::string &filename,
			      bool write_as_cif_flag,
			      mmdb::byte gz,
			      bool write_hydrogens = 1,  // optional arg
			      bool write_aniso_records = 1,  // optional arg
			      bool write_conect_records = 0);  // optional arg

class access_mol : public mmdb::Manager {
 public:
   // we use a pointer so that the destuctor doesn't
   // get run
   const mmdb::Title  *GetTitle()  const {return &title; }
};

class access_title : public mmdb::Title {
 public:
   // we use a pointer so that the destuctor doesn't
   // get run
   mmdb::TitleContainer *GetCompound() { return &compound; }
   mmdb::TitleContainer *GetAuthor() { return &author; }
};

namespace coot {

   // used by above
  void delete_hydrogens_from_mol(mmdb::Manager *mol);
  void delete_aniso_records_from_atoms(mmdb::Manager *mol);
  
  std::vector<std::string> get_compound_lines(mmdb::Manager *mol);
  std::string get_title(mmdb::Manager *mol);
}


// needs <iostream>
// 
std::ostream& operator<<(std::ostream& s, mmdb::Atom &atom);

std::ostream& operator<<(std::ostream& s, mmdb::PAtom atom); 

#endif // MMDB_H
