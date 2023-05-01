/* coords/mmdb-extras.h
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

#ifndef MMDB_EXTRAS_H
#define MMDB_EXTRAS_H

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

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

#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/contact-info.hh"

// debug this struct
void
debug_atom_selection_container(atom_selection_container_t asc);

mmdb::LinkContainer empty_links_container();


// Bond things
//
// enum bond_colours { green, red, blue, yellow, white, grey }; 
// GREEN Cl and F
// DARK_BROWN is for Br
// ORANGE is P
// DARK_GREEN is Mg
// DARK_ORANGE is Fe
enum bond_colours { CARBON_BOND, YELLOW_BOND, RED_BOND, BLUE_BOND, GREEN_BOND, MAGENTA_BOND,
		    GREY_BOND, ORANGE_BOND, CYAN_BOND, HYDROGEN_GREY_BOND,
		    DARK_BROWN_BOND, DARK_GREEN_BOND, DARK_ORANGE_BOND, DEUTERIUM_PINK,
                    DARK_VIOLET, // I
                    VIOLET // Li, NA, K, Rb, Cs, Fr
};

float max_bond_length(const std::string &element);

int get_atom_colour_from_element(const std::string &element); 

// CCP4 symmetry library checking:
int check_ccp4_symm(); 

namespace coot { 
  // 
  // Note, we also create a chain and add this residue to that chain.
  // We do this so that we have a holder for the segid.
  // 
  // if atom_index_handle is not negative, then we try to copy the
  // "atoms index" udd atom_indices to "old atom index".  If it is, we
  // don't.
  // 
  // whole_residue_flag: only copy atoms that are either in this altLoc,
  // or has an altLoc of "".
  //
  // Returns NULL if input residue has no atoms
  // 
  // caller should delete the chain of this residue (which will
  // implicitly delete this residue too).
  // 
  mmdb::Residue *
  deep_copy_this_residue_old_style(mmdb::Residue *residue,
                                   const std::string &altconf,
                                   short int whole_residue_flag,
                                   int atom_index_handle,
                                   bool embed_in_chain_flag);

  std::pair<mmdb::Residue *, atom_selection_container_t>
    deep_copy_this_residue_and_make_asc(mmdb::Manager *orig_mol,
					mmdb::Residue *residue,
					const std::string &altconf, 
					short int whole_residue_flag,
					int atom_index_handle,
					int udd_afix_handle);

  // 13 14 15 20 21 22  -> 1
  // 13 14 15 20 22 21  -> 0
  short int progressive_residues_in_chain_check(const mmdb::Chain *chain_p); 

  // Typically this is used on an asc (moving atoms) to get the N of a
  // peptide (say).  Return NULL on atom not found.
  mmdb::Atom *get_first_atom_with_atom_name(const std::string &atomname, 
				       const atom_selection_container_t &asc); 

  // tinker with asc
  void add_atom_index_udd_as_old(atom_selection_container_t asc);


}



#endif  // MMDB_EXTRAS_H
