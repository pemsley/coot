/* ligand/monomer-utils.cc
 * 
 * Copyright 2002, 2003, 2004, 2005 by The University of York
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
 * 02110-1301, USA.
 */

#include <string.h>  // for strcpy

#include "compat/coot-sysdep.h"

#include "monomer-utils.hh"


mmdb::Residue *
coot::deep_copy_residue(mmdb::Residue *residue) {

   mmdb::Residue *rres = new mmdb::Residue;
   mmdb::Chain   *chain_p = new mmdb::Chain;
   chain_p->SetChainID(residue->GetChainID());
   rres->SetResID(residue->GetResName(),
		  residue->GetSeqNum(),
		  residue->GetInsCode());

   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   residue->GetAtomTable(residue_atoms, nResidueAtoms);
   mmdb::Atom *atom_p;
   
   for(int iat=0; iat<nResidueAtoms; iat++) { 
      atom_p = new mmdb::Atom;
      atom_p->Copy(residue_atoms[iat]);
      // std::cout << "DEBUG:: " << atom_p << std::endl;
      rres->AddAtom(atom_p);
   }
   chain_p->AddResidue(rres);
   return rres;
}

			     
void
coot::monomer_utils::add_torsion_bond_by_name(const std::string &atom_name_1,
					      const std::string &atom_name_2) {

   atom_name_pair_list.push_back(coot::atom_name_pair(atom_name_1,
						      atom_name_2)); 
}

void
coot::monomer_utils::add_torsion_bond_by_name(const std::string &atom_name_1,
					      const std::string &atom_name_2,
					      const std::string &atom_name_3,
					      const std::string &atom_name_4) {

   atom_name_quad_list.push_back(coot::atom_name_quad(atom_name_1,
						      atom_name_2,
						      atom_name_3,
						      atom_name_4)); 
}



coot::contact_info
coot::monomer_utils::getcontacts(const atom_selection_container_t &asc) const {

   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   float min_dist = 0.1;
   float max_dist = 1.9; // CB->SG CYS 1.8A
   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   asc.mol->SeekContacts(asc.atom_selection, asc.n_selected_atoms,
			 asc.atom_selection, asc.n_selected_atoms,
			 min_dist, max_dist, // min, max distances
			 0,        // seqDist 0 -> in same res also
			 pscontact, n_contacts,
			 0, &my_matt, i_contact_group);

   return contact_info(pscontact, n_contacts);

} 


std::vector<coot::atom_index_pair> 
coot::monomer_utils::get_atom_index_pairs(const std::vector<coot::atom_name_pair> &atom_name_pairs_in,
				     const mmdb::PPAtom atoms, int nresatoms) const {

   int i_store_index;
   std::vector<coot::atom_index_pair> index_pairs;

   for (unsigned int ipair=0; ipair<atom_name_pairs_in.size(); ipair++) {
      int ifound = 0;
      i_store_index = -1;
      for(int i=0; i<nresatoms; i++) {
	 std::string atomname = atoms[i]->name;
	 if (atomname == atom_name_pairs_in[ipair].atom1) {
	    i_store_index = i;
	 }
      }
      if (i_store_index > -1) { // i.e. we found the first atom
	 for(int i2=0; i2<nresatoms; i2++) {
	    std::string atomname = atoms[i2]->name;
	    if (atomname == atom_name_pairs_in[ipair].atom2) {
	       index_pairs.push_back(coot::atom_index_pair(i_store_index, i2));
	    }
	 }
      } else {
	 std::cout << "first atom " << atom_name_pairs_in[ipair].atom1
		   << " not found in residue\n";
      }
   }

   if (index_pairs.size() != atom_name_pairs_in.size()) {
      std::cout << "failure to find all atom pair in residue atoms\n" ;
   } else {
      // std::cout << "DEBUG:: found all pairs in residue atoms\n" ;
   } 
   return index_pairs;
}


std::vector<coot::atom_index_quad>
coot::monomer_utils::get_atom_index_quads(const std::vector<coot::atom_name_quad> &atom_name_quads_in,
					  const mmdb::PPAtom atoms, int nresatoms) const {

   std::vector<coot::atom_index_quad> v;
   for (unsigned int iquad=0; iquad<atom_name_quads_in.size(); iquad++) {
      int ifound = 0;
      for (int i1=0; i1<nresatoms; i1++) {
	 std::string atom_name = atoms[i1]->name;
	 if (atom_name == atom_name_quads_in[iquad].atom_name(0)) {
	    for (int i2=0; i2<nresatoms; i2++) {
	       std::string atom_name = atoms[i2]->name;
	       if (atom_name == atom_name_quads_in[iquad].atom_name(1)) {
		  for (int i3=0; i3<nresatoms; i3++) {
		     std::string atom_name = atoms[i3]->name;
		     if (atom_name == atom_name_quads_in[iquad].atom_name(2)) {
			for (int i4=0; i4<nresatoms; i4++) {
			   std::string atom_name = atoms[i4]->name;
			   if (atom_name == atom_name_quads_in[iquad].atom_name(3)) {
			      v.push_back(coot::atom_index_quad(i1, i2, i3, i4));
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   if (v.size() < atom_name_quads_in.size()) {
      std::cout << "Monomer utils: Failure to find correct atom quads in residue atoms\n" ;
      for (unsigned int iquad=0; iquad<atom_name_quads_in.size(); iquad++) {
	 std::cout << "  quad needed: :"
		   << atom_name_quads_in[iquad].atom_name(0) << ":  :"
		   << atom_name_quads_in[iquad].atom_name(1) << ":  :"
		   << atom_name_quads_in[iquad].atom_name(2) << ":  :"
		   << atom_name_quads_in[iquad].atom_name(3) << ":\n";
      }
      for (unsigned int iv=0; iv<v.size(); iv++) {
	 std::cout << "  found quad: "
		   << v[iv].index1 << "  "
		   << v[iv].index2 << "  "
		   << v[iv].index3 << "  "
		   << v[iv].index4 << "\n";
      }
      for (int i1=0; i1<nresatoms; i1++) {
	 std::cout << "  res atom " << i1 << " " << atoms[i1] << "\n";
      }
   } else {
      // std::cout << "found all quads in residue atoms\n" ;
   } 
   return v;
}

// [1-indexed]
std::vector<std::pair<int, float> >
coot::monomer_utils::get_chi_angles(mmdb::Residue *residue) const {
   
   std::vector<std::pair<int, float> > v;
   std::vector<coot::atom_index_quad> quads = get_quads(atom_name_quad_list, residue);
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);
   
   for (unsigned int i_quad = 0; i_quad<quads.size(); i_quad++) {
      clipper::Coord_orth p1(atom_to_co(residue_atoms[quads[i_quad].index1]));
      clipper::Coord_orth p2(atom_to_co(residue_atoms[quads[i_quad].index2]));
      clipper::Coord_orth p3(atom_to_co(residue_atoms[quads[i_quad].index3]));
      clipper::Coord_orth p4(atom_to_co(residue_atoms[quads[i_quad].index4]));
      double tors =  clipper::Util::rad2d(clipper::Coord_orth::torsion(p1, p2, p3, p4));
      v.push_back(std::pair<int,float>(i_quad+1, tors));
   }
   return v;
}

std::vector<coot::atom_index_quad>
coot::monomer_utils::get_quads(const std::vector<coot::atom_name_quad> &atom_name_quads,
			       mmdb::Residue *residue) const {
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);
   return get_atom_index_quads(atom_name_quads, residue_atoms, n_residue_atoms);
}


// files in ligand directory cannot depend on coords headers!

// coot::Cartesian
// coot::monomer_utils::coord_orth_to_cartesian(const clipper::Coord_orth &c) {
//    return coot::Cartesian(c.x(), c.y(), c.z());
// }

// clipper::Coord_orth coot::monomer_utils::coord_orth_to_cart(const Cartesian &c) {
//    return clipper::Coord_orth(c.x(), c.y(), c.z());
// } 

clipper::Coord_orth
coot::monomer_utils::atom_to_co(mmdb::Atom *at) const {
   return clipper::Coord_orth(at->x, at->y, at->z);
}
