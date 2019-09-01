/* src/get-monomer.hh
 * 
 * Copyright 2010 by the University of Oxford
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

#include <string>

/*! \file
  \brief Coot Scripting Interface - Monomer Interface
*/

/*! \brief import libcheck monomer give the 3-letter code. 

@return the new molecule number, if not -1 (error). */
int get_monomer(const std::string &comp_id);

//! get the monomer for the given molecule
int get_monomer_for_molecule(const std::string &comp_id, int imol);

/*! \brief Use the protein geometry dictionary to retrieve a set of
   coordinates quickly from cif data read in from the RCSB's Chemical
   Component Library.  There are no restraints from this method
   though. */
int get_monomer_from_dictionary(const std::string &comp_id, int idealised_flag);


// return a new molecule number
int get_monomer_molecule_by_network_and_dict_gen(const std::string &text);

// int get_monomer_for_molecule_by_index() is called in callbacks.h, needs a c-interface.
