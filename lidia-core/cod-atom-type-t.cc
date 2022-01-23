/* lidia-core/cod-atom-type-t.cc
 * 
 * Copyright 2016 by Medical Research Council
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

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <algorithm>

#include <sys/time.h> // for time testing

#include "rdkit-interface.hh"
#include "GraphMol/Rings.h"
#include "GraphMol/RingInfo.h"
#include "boost/dynamic_bitset.hpp"

#include "utils/coot-utils.hh"

#include "primes.hh"
#include "cod-atom-types.hh"


// we are making cod types for which base atom?
// base_atom_p.
//
// neighbour info atom type: e.g. L3: "N(CCS)(SN)" needs L2: "C-3:S-2:"
//
// cod::atom_types_t::make_cod_level_2_atom_type(RDKit::Atom *base_atom_p,
// 					      const RDKit::ROMol &rdkm) {


cod::atom_level_2_type::atom_level_2_type(const RDKit::Atom *base_atom_p,
					  const RDKit::ROMol &rdkm) {

   //
   // first is int,string, number of rings and atom-ring,string
   // e.g. 2,"C[5a,6a]"  1,"N[6a]"
   // 
   // second is the  vector of hybridization of the neighbours, e.g.(2,2,3)
   //
   // 20170606 this is now a member data item
   // std::vector<atom_level_2_component_type> components;
   //
   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

   RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
   boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(base_atom_p);
   while(nbrIdx != endNbrs) {
      const RDKit::Atom *at_neighb = rdkm[*nbrIdx];
      unsigned int degree = at_neighb->getDegree();
      int n = at_neighb->getAtomicNum();
      std::pair<int,std::string> ring_info = make_ring_info_string(at_neighb);
      
      std::string atom_ele = tbl->getElementSymbol(n);

      // std::string s = atom_ele;
      // s += ring_info.second;
      // s += "-";
      
      atom_level_2_component_type c(at_neighb, rdkm);
      components.push_back(c);
      
      nbrIdx++;
   }

   if (false) { // debug
      std::string atom_name;
      base_atom_p->getProp("name", atom_name);
      std::cout << "--- about to sort components for atom " << atom_name << std::endl;
   }


   std::sort(components.begin(), components.end(), level_2_component_sorter);

   if (false) { // debug

      std::string atom_name;
      base_atom_p->getProp("name", atom_name);
      // C[5a,6a] comes before C[5a]
      // C5[a]    comes before C[6a]
      //
      std::cout << "   components sorted: " << atom_name << std::endl;
      for (unsigned int i=0; i<components.size(); i++)
	 std::cout << " component  " << i << "   " << components[i] << std::endl;
   }

   // todo:
   // N with degree 3 in aromatic ring or sp2 ring should have n_extra_elect = 2
   //
   std::string name;
   base_atom_p->getProp("name", name);
   int base_valence = tbl->getDefaultValence(base_atom_p->getAtomicNum());
   int e_valence = base_atom_p->getExplicitValence();
   n_extra_elect = e_valence - base_atom_p->getDegree() + base_atom_p->getFormalCharge();

   if (false)
      std::cout << "debug-elec " << name
		<< " base-val: " << base_valence
		<< " e_val " << e_valence
		<< " degree: " << base_atom_p->getDegree()
		<< " totalDegree: " << base_atom_p->getTotalDegree()
		<< " formal-charge: " << std::setw(2) << base_atom_p->getFormalCharge()
		<< " totalvalence " << base_atom_p->getTotalValence() << " "
		<< " ->  " << n_extra_elect << "\n"
		<< std::endl;

   std::string l2;

   int i_last_component = components.size() -1;
   int n_components = components.size();
   for (int i=0; i<n_components; i++) {
      l2 += components[i].element;
      l2 += components[i].ring_info_string;
      if (components[i].neighb_degrees.size()) {
	 l2 += "-";
	 int n_size = components[i].neighb_degrees.size(); // int for comparison below
	 for (int j=0; j<n_size; j++) {
	    if (j != 0)
	       l2 += "_";
	    l2 += coot::util::int_to_string(components[i].neighb_degrees[j]);
	 }
	 if (i != i_last_component)
	    l2 += ":";
      }
   }

   str = l2;
}

std::string
cod::atom_level_2_type::extra_electron_type() const {

   std::string s;

   int i_last_component = components.size() -1;
   int n_components = components.size();
   for (int i=0; i<n_components; i++) {
      s += components[i].element;
      // std::cout << "s: " << s << std::endl;
      s += components[i].ring_info_string;
      // std::cout << "s: " << s << std::endl;
      if (components[i].neighb_degrees.size()) {
	 s += "-";
	 // std::cout << "s: " << s << std::endl;
	 int n_size_1 = components[i].neighb_degrees.size(); // int for comparison below
	 int n_size = components[i].neighb_extra_elect.size(); // int for comparison below
	 if (n_size != n_size_1) {
	    std::cout << "------ mismatch sizes: " << n_size << " " << n_size_1 << std::endl;
	 }
	 int idx_last = n_size - 1;
	 for (int j=0; j<n_size; j++) { 
	    if (j != 0) {
	       s += "_";
	       // std::cout << "s: " << s << std::endl;
	    }
	    int nee = components[i].neighb_extra_elect[j];
	    s += coot::util::int_to_string(nee);
	 }
	 if (i != i_last_component) {
	    s += ":";
	    // std::cout << "s: " << s << std::endl;
	 }
      }
   }
   return s;
}

int
cod::atom_level_2_type::n_extra_electrons() const {

   return n_extra_elect;
}

void
cod::atom_type_t::set_neighb_degrees_string() {

   std::string s;
   for (unsigned int i=0; i<neighb_degrees.size(); i++) {
      s += coot::util::int_to_string(neighb_degrees[i]);
      s += ":";
   }
   neighb_degrees_str_ = s;
}


std::ostream &
cod::operator<<(std::ostream &s,
		const atom_level_2_type::atom_level_2_component_type &c) {

   s << "{" << c.element << " n-rings: " << c.number_of_rings << " ";
   if (! c.ring_info_string.empty())
       s << "ring-info: " << c.ring_info_string << " ";
   // neighb_degrees is a vector of ints
   if (! c.neighb_degrees.empty())
      s << "neighb-degrees ";
   for (unsigned int i=0; i<c.neighb_degrees.size(); i++)
      s << c.neighb_degrees[i] << " ";
   s << "}";
   return s;

}



//
// maybe this is the wrong file for this function
//
// static
std::string
cod::atom_type_t::level_4_type_to_level_3_type(const std::string &l4t)  {

   std::string s = l4t;

   if (! s.empty()) {

      // s.back() is C++-11
      //
      if (*(s.rbegin()) == '}') {
         // trim off {xxxxx} from s
         std::string::size_type ii = s.find_last_of('{');
         if (ii != std::string::npos)
            s = s.substr(0,ii);
      }
   }
   return s;
}


//
cod::atom_type_t::atom_type_t(const std::string &s1_hash,
			      const std::string &colon_degrees_type,
			      const cod::atom_level_2_type &l2,
			      const std::string &s3, const std::string &s4) {
   level_2 = l2;
   level_3 = s3;
   level_4 = s4;
   neighb_degrees_str_ = colon_degrees_type;

   try {
      hash_value = coot::util::string_to_int(s1_hash);
   }
   catch (const std::runtime_error &rte) {
      hash_value = -1;
   }
}


int
cod::hybridization_to_int(RDKit::Atom::HybridizationType h) {

   int hybrid = 0;
   if (h == RDKit::Atom::SP)
      hybrid = 1;
   if (h == RDKit::Atom::SP2)
      hybrid = 2;
   if (h == RDKit::Atom::SP3)
      hybrid = 3;

   return hybrid;
}

// static
bool
cod::atom_level_2_type::level_2_component_sorter(const atom_level_2_component_type &la,
						 const atom_level_2_component_type &lb) {

   // sp3 before sp2
   //
   // in 0 rings before in 1 or more rings
   //
   // in 2 rings before in 1 ring
   //
   // 5-mem ring before 6-mem ring
   //
   // neighb hybrids: _2_0 before _2_2
   //
   // C[5a]-2_1_1 before C[5a]-1_0_0 

   // std::cout << "comparing l2 components: " << la << " " << lb << std::endl;

   auto are_different_vectors = [] (const std::vector<int> &nd_1, const std::vector<int> &nd_2) {
                                   if (nd_1.size() != nd_2.size()) return false;
                                   if (nd_1.empty()) return true;
                                   bool all_same = true;
                                   for (unsigned int i=0; i<nd_1.size(); i++) {
                                      if (nd_1[i] != nd_2[i]) {
                                         all_same = false;
                                         break;
                                      }
                                   }
                                   return ! all_same;
                                };

   bool extra_electron_sort = true; // 20170613: inner sort by extra electrons, not hybridization

   bool status = false;
   if (la.element < lb.element) {
      return true;
   } else {
      if (la.element > lb.element) {
	 return false;
      } else {

	 // same element

	       if ((la.number_of_rings == 0) && (lb.number_of_rings > 0)) {
		  return true;
	       } else {
		  if ((la.number_of_rings > 0) && (lb.number_of_rings == 0)) {
		     return false;
		  } else { 

		     if (la.number_of_rings > lb.number_of_rings) {
			return true;
		     } else {
			if (la.number_of_rings < lb.number_of_rings) {
			   return false;
			} else {

			   // same number of rings

			   // std::cout << "comparing ring info strings: "
			   // << la.ring_info_string << " " << lb.ring_info_string << std::endl;

			   // aromatic rings come first, "6a" is longer than "6"
			   if (lb.ring_info_string.length() > la.ring_info_string.length()) {
			      return false;
			   } else {
			      if (lb.ring_info_string.length() < la.ring_info_string.length()) {
				 return true;
			      } else {

				 if (la.ring_info_string < lb.ring_info_string) {
				    return true;
				 } else {
				    if (la.ring_info_string > lb.ring_info_string) {
				       return false;
				    } else {

				       bool s = (lb.neighb_degrees.size() <
						 la.neighb_degrees.size());

				       if (s) {
					  return true;
				       } else {

					  bool sr = (lb.neighb_degrees.size() >
						     la.neighb_degrees.size());

					  if (sr) {
					     return false;
					  } else {

                                             // so la.neighb_degrees is the same length as lb.neighb_degrees

                                             if (are_different_vectors(la.neighb_degrees, lb.neighb_degrees)) {

                                                const auto &nd_1 = la.neighb_degrees;
                                                const auto &nd_2 = lb.neighb_degrees;

                                                for (unsigned int i=0; i<nd_1.size(); i++) {
                                                   if (nd_1[i] > nd_2[i]) return true;
                                                   if (nd_1[i] < nd_2[i]) return false;
                                                }

                                             } else {

                                                if (extra_electron_sort) {

                                                   unsigned int n = la.neighb_degrees.size();
                                                   for (unsigned int i=0; i<n; i++) {
                                                      if (la.neighb_extra_elect[i] < lb.neighb_extra_elect[i]) {
                                                         return true;
                                                      } else {
                                                         if (la.neighb_extra_elect[i] > lb.neighb_extra_elect[i]) {
                                                            return false;
                                                         }
                                                      }
                                                   }
                                                   return false; // or something.

                                                } else {

                                                   // old style hybridizations

                                                   // so la.neighb_degrees should be the same
                                                   // size as lb.neighb_degrees.

                                                   unsigned int n_h_a = la.neighb_degrees.size();
                                                   for (unsigned int i=0; i<n_h_a; i++) {
                                                      if (la.neighb_degrees[i] < lb.neighb_degrees[i])
                                                         return true;
                                                      if (la.neighb_degrees[i] > lb.neighb_degrees[i])
                                                         return false;
                                                   }

                                                   return false; // or something.
                                                }
                                             }
                                          }
				       }
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
      }
   }
   return status;
}


// return number_of_ring,atom_ring_string
std::pair<int, std::string>
cod::make_ring_info_string(const RDKit::Atom *atom_p) {

   std::string atom_ring_string;
   
   std::vector<int> ring_size_vec;
   std::vector<int> ring_is_arom_vec; // flags

   try {
      atom_p->getProp("ring_size", ring_size_vec);
      atom_p->getProp("ring_arom", ring_is_arom_vec);

      // sort ring_info so that the rings with more atoms are at
      // the top.  Practically 6-rings should come above 5-rings.

      // sorting by ring size is done when the ring info is constructed
      // and added to atoms - no here
      // std::sort(ring_si.begin(), ring_info.end());
	 
      atom_ring_string = "[";
      for (unsigned int i_ring=0; i_ring<ring_size_vec.size(); i_ring++) {

	 if (i_ring > 0)
	    atom_ring_string += ",";

	 // If this ring is aromatic we need at append "a" to the (typically)
	 // n_atoms_in_ring = 6.
	 // 
	 int n_atoms_in_ring = ring_size_vec[i_ring];

	 atom_ring_string += coot::util::int_to_string(n_atoms_in_ring);
	 int arom = ring_is_arom_vec[i_ring];
	 if (arom) {
	    atom_ring_string += "a";
	 }
      }
      atom_ring_string += "]";
   }
   catch (const KeyErrorException &kee) {
      // no ring info on that atom, that's fine
      if (false)
	 std::cout << "  debug-info:: no ring info for atom " << atom_p
		   << std::endl;
   }
   return std::pair<int,std::string> (ring_size_vec.size(), atom_ring_string);
}


cod::atom_level_2_type::atom_level_2_component_type::atom_level_2_component_type(const RDKit::Atom *at, const RDKit::ROMol &rdkm) {
      
   // Version 147 atom types
   // s += coot::util::int_to_string(degree);

   // 152+ we also need info on the hybridization of the neighbours of at_neighb
   //
   std::vector<RDKit::Atom::HybridizationType> hv;
   RDKit::ROMol::ADJ_ITER nbrIdx_n, endNbrs_n;
   boost::tie(nbrIdx_n, endNbrs_n) = rdkm.getAtomNeighbors(at);
   while(nbrIdx_n != endNbrs_n) {
      const RDKit::Atom *at_neighb_n = rdkm[*nbrIdx_n];
      hv.push_back(at_neighb_n->getHybridization());
      int n_extra_elect = at_neighb_n->getExplicitValence() - at_neighb_n->getDegree() + at_neighb_n->getFormalCharge();
      neighb_extra_elect.push_back(n_extra_elect);
      nbrIdx_n++;
   }
   std::pair<int,std::string> ring_info = make_ring_info_string(at);

   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   int base_valence = tbl->getDefaultValence(at->getAtomicNum());
   int e_valence = at->getExplicitValence();
   n_extra_elect = e_valence - at->getDegree() + at->getFormalCharge();

   number_of_rings = ring_info.first;
   ring_info_string = ring_info.second;

   std::vector<int> v(hv.size());
   for (unsigned int jj=0; jj<hv.size(); jj++)
      v[jj] = hybridization_to_int(hv[jj]);
   std::sort(v.begin(), v.end());
   std::reverse(v.begin(), v.end());

   // std::cout << "v size: " << v.size() << std::endl;
   // for (unsigned int ii=0; ii<v.size(); ii++)
   // std::cout << "post sort " << v[ii] << std::endl;

   // std::pair<std::pair<int, std::string>, std::vector<int> > c(ring_info,v);

   element = tbl->getElementSymbol(at->getAtomicNum());

   neighb_degrees = v;

   // neighb extra electrons, decreasing order
   std::sort(neighb_extra_elect.begin(),
	     neighb_extra_elect.end());
   std::reverse(neighb_extra_elect.begin(),
		neighb_extra_elect.end());

   // debugging, let's also extract the atom name:
   if (true) {
      at->getProp("name", atom_name);
   }

}


#endif // MAKE_ENHANCED_LIGAND_TOOLS
