/* coords/atom-selection-container.hh
 * -*-c++-*-  
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

// mmdb-extras

// needs mmdb_manager.h and <string>

#ifndef ATOM_SELECTION_CONTAINER_HH
#define ATOM_SELECTION_CONTAINER_HH

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <iostream>

#include "math.h"

#ifndef MMDB_MANAGER_H
#define MMDB_MANAGER_H
#include <mmdb2/mmdb_manager.h>
#endif

#include "geometry/protein-geometry.hh"

// we need this for bonded_pair_container_t.
#include "coot-utils/bonded-pairs.hh"
// and atom quads
#include "mini-mol/atom-quads.hh"

// 
//
class atom_selection_container_t {

   void fill_links(mmdb::Manager *mol_other);
   void fill_links() { fill_links(mol); } 
   
public:
   mmdb::Manager *mol;
   int n_selected_atoms; 
   mmdb::PPAtom atom_selection; 
   std::string read_error_message;
   int read_success;
   int SelectionHandle;
   int UDDAtomIndexHandle; // no negative is OK.
   int UDDOldAtomIndexHandle; // ditto. // set initially to -1 in make_asc

   std::vector<mmdb::Link> links;

   atom_selection_container_t(mmdb::Manager *mol_in, int selhnd) : atom_selection(0) {

     mol = mol_in;
     mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);
     SelectionHandle = selhnd;
     UDDAtomIndexHandle = -1;
     UDDOldAtomIndexHandle = -1;
     read_success = 1;
     fill_links();
   }

   atom_selection_container_t(mmdb::PPAtom residue_atoms, int n_residue_atoms) : atom_selection(residue_atoms) {
      mol = NULL;;
      n_selected_atoms = n_residue_atoms;
      SelectionHandle = -1;
      UDDAtomIndexHandle = -1;
      UDDOldAtomIndexHandle = -1;
      read_success = 1;
   }

   atom_selection_container_t() : mol(0), n_selected_atoms(0), atom_selection(0) {
      mol = NULL;
      SelectionHandle = -1;
      UDDAtomIndexHandle = -1;
      UDDOldAtomIndexHandle = -1;
      read_success = 0;
   }

   bool empty() const {
      return (mol == NULL);
   }

   void fill_links_using_mol(mmdb::Manager *mol_other) {
      fill_links(mol_other);
   } 
  
   void clear_up() { 
      if (read_success) 
	 if (SelectionHandle)
	    mol->DeleteSelection(SelectionHandle);
      delete mol;
      atom_selection = 0;
      mol = 0;
   }

   void delete_atom_selection() { 
      mol->DeleteSelection(SelectionHandle);
      n_selected_atoms = 0;
      atom_selection = 0;
   }

   // sets UDDOldAtomIndexHandle
   void add_old_atom_indices();

   //! return 0,0,0 on no centre.
   //! centre is not weighted by occupancy or atomic number
   clipper::Coord_orth get_selection_centre() const {
      clipper::Coord_orth sum(0,0,0);
      unsigned int count = 0;
      for (int i=0; i<n_selected_atoms; i++) {
	 mmdb::Atom *at = atom_selection[i];
	 if (at) {
	    sum += clipper::Coord_orth(at->x, at->y, at->z);
	    count++;
	 }
      }
      if (count > 0) {
	 float d = 1.0/static_cast<float>(count);
	 sum = d * sum;
	 return sum;
      } else {
	 return sum;
      }
   }
   
   void apply_shift(float x_shift, float y_shift, float z_shift) {
      for (int i=0; i<n_selected_atoms; i++) { 
	 atom_selection[i]->x -= x_shift;
	 atom_selection[i]->y -= y_shift;
	 atom_selection[i]->z -= z_shift;
      }
   }

   mmdb::Residue *get_next(mmdb::Residue *) const;
   mmdb::Residue *get_previous(mmdb::Residue *) const;
};

atom_selection_container_t make_asc(mmdb::Manager *mol, bool transfer_atom_indices_flag=false);

atom_selection_container_t get_atom_selection(std::string t, 
					      bool allow_duplseqnum,
                                              bool verbose_mode,
					      bool convert_to_v2_name_flag);

// put these in coot namespace? -- FIXME

void fix_element_name_lengths(mmdb::Manager *mol);
int fix_nucleic_acid_residue_names(atom_selection_container_t asc);
int fix_nucleic_acid_residue_name(mmdb::Residue *r); // return whether it was changed or not.
void convert_to_old_nucleotide_atom_names(mmdb::Residue *r);

// return the number of fixed atoms
int fix_away_atoms(atom_selection_container_t asc);

// return the number of changed hydrogen names
int fix_wrapped_names(atom_selection_container_t asc);

// OK, we don't wan't lidia-core functions in coot-utis
// (lidia-core depends on coot-utils)
//
// #include "lidia-core/lbg-molfile.hh"
// #include "lidia-core/lig-build.hh"

namespace coot { 
   bool is_hydrogen(const std::string &ele);
   bool is_deuterium(const std::string &ele);
}

void debug_atom_selection_container(atom_selection_container_t asc);


#endif // ATOM_SELECTION_CONTAINER_HH
