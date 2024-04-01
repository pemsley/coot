/*
 * geometry/protein-donor-acceptors.cc
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

#include <iostream>

#include "protein-donor-acceptors.hh"

void
coot::quick_protein_donor_acceptors::init() {

   std::pair<std::string, std::string> p;

   // PDBv3 FIXME
   hb_type_map[key("CYS", " SG ")] = HB_ACCEPTOR;
   hb_type_map[key("ASP", " OD1")] = HB_ACCEPTOR;
   hb_type_map[key("ASP", " OD2")] = HB_ACCEPTOR;
   hb_type_map[key("GLU", " OE1")] = HB_ACCEPTOR;
   hb_type_map[key("GLU", " OE2")] = HB_ACCEPTOR;
   hb_type_map[key("HIS", " ND1")] = HB_BOTH;
   hb_type_map[key("HIS", " NE2")] = HB_BOTH;
   hb_type_map[key("LYS", " NZ ")] = HB_DONOR;
   hb_type_map[key("MET", " SD ")] = HB_ACCEPTOR;
   hb_type_map[key("ASN", " OD1")] = HB_ACCEPTOR;
   hb_type_map[key("ASN", " ND2")] = HB_DONOR;
   hb_type_map[key("GLN", " OE1")] = HB_ACCEPTOR;
   hb_type_map[key("GLN", " NE2")] = HB_DONOR;
   hb_type_map[key("ARG", " NE ")] = HB_DONOR;
   hb_type_map[key("ARG", " NH1")] = HB_DONOR;
   hb_type_map[key("ARG", " NH2")] = HB_DONOR;
   hb_type_map[key("SER", " OG ")] = HB_BOTH;
   hb_type_map[key("THR", " OG ")] = HB_BOTH;
   hb_type_map[key("TRP", " NE1")] = HB_DONOR;
   hb_type_map[key("TYR", " OH ")] = HB_BOTH;

   const char *l[] = {"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
		      "MET", "MSE", "ASN", "PRO", "GLN", "ARG", "SER", "TYR", "VAL", "TRP",
		      "TYR" };
   unsigned int n_res_types = 21;

   for (unsigned int i=0; i<n_res_types; i++) {
      key pn(l[i], " N  ");
      key po(l[i], " O  ");
      hb_type_map[pn] = HB_DONOR;
      hb_type_map[po] = HB_ACCEPTOR;
   }

}

coot::hb_t
coot::quick_protein_donor_acceptors::get_type(const key &k) const {

   hb_t r = HB_UNASSIGNED;

   std::map<key, hb_t>::const_iterator it = hb_type_map.find(k);

   if (it != hb_type_map.end()) {
      r = it->second;
   }

   return r;
}

// first: did we find an answer (ie. both were protein residues?)
// second: was it a potential hydrogen bond (had the correct atom types?"
//
std::pair<bool, bool>
coot::quick_protein_donor_acceptors::is_hydrogen_bond_by_types(const std::pair<key, key> &hbtp) const {

   bool found = false;
   bool r = false;
   const key &key_1 = hbtp.first;
   const key &key_2 = hbtp.second;
   hb_t type_1 = get_type(key_1);
   if (type_1 == HB_BOTH || type_1 == HB_ACCEPTOR || type_1 == HB_DONOR) {
      hb_t type_2 = get_type(key_2);
      if (type_2 == HB_BOTH || type_2 == HB_ACCEPTOR || type_2 == HB_DONOR) {
	 found = true;
	 if (type_1 == HB_BOTH || type_1 == HB_ACCEPTOR)
	    if (type_2 == HB_BOTH || type_2 == HB_DONOR)
	       r = true;
	 if (type_1 == HB_BOTH || type_1 == HB_DONOR)
	    if (type_2 == HB_BOTH || type_2 == HB_ACCEPTOR)
	       r = true;
      } else {
	 if (type_2 == HB_NEITHER)
	    found = true;
      }
   } else {
      hb_t type_2 = get_type(key_2);
      if (type_1 != HB_UNASSIGNED)
	 if (type_2 != HB_UNASSIGNED)
	    found = true;
   }

   return std::pair<bool, bool> (found, r);
}

std::vector<std::pair<bool, bool> >
coot::quick_protein_donor_acceptors::is_hydrogen_bond_by_types(std::vector<std::pair<key, key> > &hbtp) const {

   std::vector<std::pair<bool, bool> > v(hbtp.size());
   for (unsigned int i=0; i<hbtp.size(); i++)
      v[i] = is_hydrogen_bond_by_types(hbtp[i]);
   return v;
}



void
coot::quick_protein_donor_acceptors::test() const {
   
   hb_t r;

   if (false) {
      std::cout << "ref: HB_UNASSIGNED " << HB_UNASSIGNED << std::endl;
      std::cout << "ref: HB_NEITHER " << HB_NEITHER << std::endl;
      std::cout << "ref: HB_DONOR " << HB_DONOR << std::endl;
      std::cout << "ref: HB_ACCEPTOR " << HB_ACCEPTOR << std::endl;
      std::cout << "ref: HB_HYDROGEN " << HB_HYDROGEN << std::endl;
   }

   bool r1 = (is_hydrogen_bond_by_types(std::pair<key,key>(key("ALA", " CB "), key("SER", " OG "))).second == false);
   bool r2 = (is_hydrogen_bond_by_types(std::pair<key,key>(key("ALA", " N  "), key("SER", " OG "))).second == true);
   bool r3 = (is_hydrogen_bond_by_types(std::pair<key,key>(key("TYR", " N  "), key("TRP", " NE1"))).second == false);

   std::cout << r1 << " " << r2 << " " << r3 << std::endl;
}
