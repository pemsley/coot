/*
 * geometry/main-chain.cc
 *
 * Copyright 2016 by Medical Research Council
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

#include "main-chain.hh"


// PDBv3 FIXME

bool
coot::is_main_chain_p(mmdb::Atom *at) { 

   std::string mol_atom_name(at->name);
   if (mol_atom_name == " N  " ||
       mol_atom_name == " C  " ||
       mol_atom_name == " CA " ||
       mol_atom_name == " H  " ||
       mol_atom_name == " HA " || // CA hydrogen
       mol_atom_name == " OXT" ||
       mol_atom_name == " O  ") {
      return true;
   } else {
      std::string res_name = at->GetResName();
      if (res_name == "GLY") {
	 if (mol_atom_name == " HA2" ||
	     mol_atom_name == " HA3") {
	    return 1;
	 }
	 
	 // Perhaps N-terminal H atom?
	 mmdb::Residue *res = at->residue;
	 if (res) {
	    if (res->isNTerminus()) {
	       if (mol_atom_name == " H1 ") return true;
	       if (mol_atom_name == " H2 ") return true;
	       if (mol_atom_name == " H3 ") return true;
	    }
	 }
      }
      return 0;
   } 
}

bool
coot::is_main_chain_or_cb_p(mmdb::Atom *at) { 

   std::string mol_atom_name(at->name);
   return is_main_chain_or_cb_p(mol_atom_name);
}

// return 0 or 1
bool
coot::is_main_chain_p(const std::string &mol_atom_name) {

   if (mol_atom_name == " N  " ||
       mol_atom_name == " C  " ||
       mol_atom_name == " H  " ||
       mol_atom_name == " CA " ||
       mol_atom_name == " CB " ||
       mol_atom_name == " HA " || // CA hydrogen
       mol_atom_name == " O  ") {
      return 1;
   } else {
      return 0;
   } 
}

// return 0 or 1
bool
coot::is_main_chain_or_cb_p(const std::string &mol_atom_name) {

   if (mol_atom_name == " N  " ||
       mol_atom_name == " C  " ||
       mol_atom_name == " H  " ||
       mol_atom_name == " CA " ||
       mol_atom_name == " OXT" ||
       mol_atom_name == " CB " ||
       mol_atom_name == " HA " || // CA hydrogen
       mol_atom_name == " O  ") {
      return 1;
   } else {
      return 0;
   } 
} 


// return 0 or 1
bool coot::is_hydrogen_p(mmdb::Atom *at) {

   std::string mol_atom_ele(at->element);
   if (mol_atom_ele == " H" ||
       mol_atom_ele == " D") {
      return 1;
   } else {
      return 0;
   }
}
