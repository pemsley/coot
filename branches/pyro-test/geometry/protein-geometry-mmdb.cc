/* geometry/protein-geometry.cc
 * 
 * Copyright 2004 by Eugene Krissinel
 * Copyright 2004, 2005 The University of York
 * 
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

#include <iostream>
#include "protein-geometry.hh"

#include <mmdb/mmdb_manager.h> // for debugging.

// This function is based on mmdb's CATom::MakePDBAtomName(), so we share
// the copyright.
// 
std::string
coot::protein_geometry::comp_atom_pad_atom_name(const std::string &atom_id,
						const std::string &type_symbol) const {

   std::string element = type_symbol;
   std::string name    = atom_id;
   std::string new_name; // = atom_id;


   if (atom_id.length() == 4) {

      new_name = atom_id;

   } else {
   
      if (element == "") {
	 if (name.length() == 1) {
	    new_name = " ";
	    new_name += name;
	    new_name += "  ";
	 }
      
      } else { // element was defined

	 // consider the following cases of 
	 // atom_id and element:
	 //   "CA"  "C"  -> " CA "
	 //   "AC8" "C"  -> "AC8 " : in this case we post apply
	 //                          the space because the second char of
	 //                          the atom_id is equal to the
	 //                          element
	 //   "NN1" "N"  -> "NN1 " :
	 //   "NP"  "P"  -> "NP  "
	 
	 if (element.length() == 1) {
	    int k = atom_id.length();
	    if (k == 3) {
	       // std::cout << "comparing " << atom_id.substr(1,1) << " and " << element << std::endl;
	       if (atom_id.substr(0,1) == element) {
		  new_name = " " + atom_id;
	       } else { 
		  if (atom_id.substr(1,1) == element) {
		     new_name = atom_id + " ";
		  } else {
		     new_name = " " + atom_id;
		  }
	       }
	    } else { 
	       // promote the characters one space
	       if (k==2) {

		  // e.g "NP" "P", or "CA" "C" (20100929, eh?)
		  if ((atom_id.substr(1,1) == element) &&
		      (element != "H") && (element != "C") && (element != "N")) {
		     new_name = atom_id + "  ";
		  } else {
		     new_name = " " + atom_id;
		     new_name += " ";
		  }
	       } else {
		  new_name = " ";
		  new_name += atom_id;
		  // fill the rest with spaces
		  new_name += "  ";
	       }
	    }

	 } else {
	    // element was 2 (or more) characters
	    if (element[0] == ' ') { // unusual from dict
	       if (element[1] != name[1]) {
		  // promote the characters one space
		  new_name = " ";
		  new_name += atom_id;
		  int k = atom_id.length();
		  // fill the rest with spaces
		  if (k == 1) 
		     new_name += "  ";
		  if (k == 2)
		     new_name += " ";
	       }
	    } else {
	       // This is the usual (always?) case, there is no leading space
	       //
	       // left justify the name and pad with spaces
	       new_name = atom_id;
	       if (atom_id.size() == 1) // never happens?
		  new_name += "   ";
	       if (atom_id.size() == 2)
		  new_name += "  ";
	       if (atom_id.size() == 3)
		  new_name += " ";
	       
	    }
	 }
      }
   }

// debug
//    
//    CAtom at;
//    at.SetAtomName(atom_id.c_str());
//    at.SetElementName(type_symbol.c_str());
//    at.MakePDBAtomName();
//    if (new_name != std::string(at.name)) {
//       std::cout << "name pad failure, mmdb coot :" << at.name << ": :"
// 		<< new_name << ": for in_atom :" << atom_id.c_str()
// 		<< ": element :" << type_symbol.c_str() << ":" << std::endl;
//    } else {
//       std::cout << "name pad match " << new_name << std::endl;
//    }

//    std::cout << "new atom name :" << new_name << ": from :"
// 	     << atom_id << ": :" << element << ":" << std::endl;
   
//    std::cout << "DEBUG:: comp_atom_pad_atom_name() for :" << atom_id
//           << ": element is :" << element
//  	     << ": returning :" << new_name << ":\n";
   
   return new_name;
}
