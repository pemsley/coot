/*
 * coords/Bond_lines_ext.h
 * 
 * Copyright 2003 by the University of York
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
 */

/* #include <string> */
/* #include <fstream> */
/* #include <vector> */

/* #include "Cartesian.h" */
/* #include "mmdb_manager.h" */
/* #include "mmdb-extras.h" */
/* #include "mmdb.h" */
/* #include "mmdb-crystal.h"  // should be merged with extras */


// #include "Bond_lines.h"

class Bond_lines_ext : public Bond_lines_container {

   // Add clipper facilities here with which we don't want to complicate
   // Bond_lines
 public:

   Bond_lines_ext() {}; // for use with find_molecule_middle(); 

   void find_skel_atom_bonds(atom_selection_container_t SelAtom);

   Bond_lines_ext(atom_selection_container_t SelAtom) {
      find_skel_atom_bonds(SelAtom); }

   coot::Cartesian find_molecule_middle(atom_selection_container_t SelAtom,
					float max_neighbour_dist); 

};
