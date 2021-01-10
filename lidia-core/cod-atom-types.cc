/* lidia-core/cod-atom-types.cc
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


void
cod::fill_element_period_group_map() {

   if (element_period_group_map.size() == 0) {
      element_period_group_map[ "H"] = std::pair<unsigned int, unsigned int>(1, 1);
      element_period_group_map[ "B"] = std::pair<unsigned int, unsigned int>(2,13);
      element_period_group_map[ "C"] = std::pair<unsigned int, unsigned int>(2,14);
      element_period_group_map[ "N"] = std::pair<unsigned int, unsigned int>(2,15);
      element_period_group_map[ "O"] = std::pair<unsigned int, unsigned int>(2,16);
      element_period_group_map[ "F"] = std::pair<unsigned int, unsigned int>(2,17);
      element_period_group_map["Si"] = std::pair<unsigned int, unsigned int>(3,14);
      element_period_group_map[ "P"] = std::pair<unsigned int, unsigned int>(3,15);
      element_period_group_map[ "S"] = std::pair<unsigned int, unsigned int>(3,16);
      element_period_group_map["Cl"] = std::pair<unsigned int, unsigned int>(3,17);
      element_period_group_map["Ge"] = std::pair<unsigned int, unsigned int>(4,14);
      element_period_group_map["As"] = std::pair<unsigned int, unsigned int>(4,15);
      element_period_group_map["Se"] = std::pair<unsigned int, unsigned int>(4,16);
      element_period_group_map["Br"] = std::pair<unsigned int, unsigned int>(4,17);
      element_period_group_map["Sb"] = std::pair<unsigned int, unsigned int>(5,15);
      element_period_group_map["Te"] = std::pair<unsigned int, unsigned int>(5,16);
      element_period_group_map[ "I"] = std::pair<unsigned int, unsigned int>(5,17);
   }
}


// can throw a std::runtime_error
//
// rdkit_mol is not const because there is no const beginAtoms() operator.
//
// add_name_as_property is an optional argument default true
// 
std::vector<cod::atom_type_t>
cod::atom_types_t::get_cod_atom_types(RDKit::ROMol &rdkm,
				      bool add_name_as_property) {


   fill_element_period_group_map();

   std::vector<atom_type_t> v;

   RDKit::RingInfo* ring_info_p = rdkm.getRingInfo();

   unsigned int n_rings = ring_info_p->numRings();

   // Maybe add a vector of ring sizes to the atoms (most atoms will
   // not have a vector added (because they are not part of rings)).
   // 
   std::vector<std::vector<int> > atomRings = ring_info_p->atomRings();

   // Now sort ring_info so that the rings with more atoms are at
   // the top.  Practically 6-rings should come after 5-rings.
   //
   std::vector<std::vector<int> > sorted_atomRings = atomRings;
   // 
   //
   std::sort(sorted_atomRings.begin(), sorted_atomRings.end(), atomRingSorter);
   
   for (unsigned int i_ring=0; i_ring<n_rings; i_ring++) {
      std::vector<int> ring_atom_indices = sorted_atomRings[i_ring];

      int arom_flag = coot::is_aromatic_ring(ring_atom_indices, rdkm); // 0 or 1
      unsigned int n_ring_atoms = ring_atom_indices.size();

      // don't include macrocycle ring info (and the like (like cycloheptane))
      if (n_ring_atoms <= 6) {
	 for (unsigned int iat=0; iat<n_ring_atoms; iat++) { 
	    try {
	       // fill these by reference
	       std::vector<int> ring_size_vec;
	       std::vector<int> ring_arom_vec;

	       // when we add a ring to ring_size, we add flag for the
	       // aromaticity also.

	       rdkm[ring_atom_indices[iat]]->getProp("ring_size", ring_size_vec);
	       ring_size_vec.push_back(n_ring_atoms);
	       rdkm[ring_atom_indices[iat]]->setProp("ring_size", ring_size_vec);

	       rdkm[ring_atom_indices[iat]]->getProp("ring_arom", ring_arom_vec);
	       ring_arom_vec.push_back(arom_flag);
	       rdkm[ring_atom_indices[iat]]->setProp("ring_arom", ring_arom_vec);
	    }
	    catch (const KeyErrorException &err) {
	       // new ring info for this atom
	       std::vector<int> ring_size_vec(1);
	       std::vector<int> ring_arom_vec(1);
	       ring_size_vec[0] = n_ring_atoms;
	       ring_arom_vec[0] = arom_flag;
	       rdkm[ring_atom_indices[iat]]->setProp("ring_size", ring_size_vec);
	       rdkm[ring_atom_indices[iat]]->setProp("ring_arom", ring_arom_vec);
	    }
	 }
      }
   }

   boost::dynamic_bitset<> fusDone(n_rings);
   unsigned int curr = 0;
   RDKit::INT_INT_VECT_MAP neighMap; // std::map<int, std::vector<int> >
   RingUtils::makeRingNeighborMap(sorted_atomRings, neighMap);
   RDKit::INT_VECT fused;
   while (curr < n_rings) {
      fused.clear();
      RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
      std::vector<std::vector<int> > fused_rings;
      fused_rings.reserve(fused.size());
      for (RDKit::INT_VECT_CI rid = fused.begin(); rid != fused.end(); ++rid) {
	 fused_rings.push_back(sorted_atomRings[*rid]);
      }

      if (fused_rings.size() > 1)
	 handle_bigger_rings_from_fused_rings(rdkm, fused_rings);
      curr++;
   }


   cod::primes primes(600000);
   
   timeval start_time;
   timeval current_time;
   // gettimeofday(&start_time, NULL);
   
   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm.beginAtoms(); ai!=rdkm.endAtoms(); ai++) {
      
      // std::pair<std::string, std::list<third_neighbour_info_t> > s_pair =
      // get_cod_atom_type(0, *ai, *ai, rdkm); // full-spec
      // std::string s = s_pair.first;
      
      atom_type_t atom_type = get_cod_atom_type(*ai, 0, 0, *ai, rdkm);

      // is this the right place?
      // 
      atom_type.set_hash_value(make_hash_index(*ai, primes));

      if (false)
	 std::cout << "atom type "
		   << atom_type.hash_value << " types: " 
		   << atom_type.level_2.string() << " "
		   << atom_type.level_4 << std::endl;
      
      v.push_back(atom_type);
      
      if (add_name_as_property)
	 (*ai)->setProp("CODAtomName", atom_type.level_4);
   }
   if (false) { // it's worth pre-computing the primes -600ms
      gettimeofday(&current_time, NULL);
      double td = current_time.tv_sec - start_time.tv_sec;  // milliseconds
      td *= 1000.0;
      td += double(current_time.tv_usec - start_time.tv_usec)/1000.0;
      std::cout << "time-diff (ms): " << td << std::endl;
   }
   
   if (v.size() != rdkm.getNumAtoms())
      throw std::runtime_error("mismatch size in get_cod_atom_types()");
      
   return v;
}

// static
bool
cod::atom_types_t::atomRingSorter(const std::vector<int> &r1, const std::vector<int> &r2) {

   return (r1.size() < r2.size());
}


void
cod::atom_types_t::handle_bigger_rings_from_fused_rings(RDKit::ROMol &rdkm,
							const std::vector<std::vector<int> > &fused_rings) {


   // first make a connection "tree".  What is bonded to what (in the
   // rings).  bond map is a list of atom indices that are bonded to
   // the key (an atom index)
   //
   std::map<int, std::vector<int> > bond_map;
   std::map<int, std::vector<int> >::iterator it;

   for (unsigned int iat=0; iat<rdkm.getNumAtoms(); iat++) {
      if (is_ring_member(iat, fused_rings)) {
	 RDKit::Atom *this_at = rdkm[iat];
	 unsigned int idx_c = this_at->getIdx();
	 RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
	 boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(this_at);
	 while(nbrIdx != endNbrs) {
	    if (is_ring_member(*nbrIdx, fused_rings)) { 
	       RDKit::Atom *at = rdkm[*nbrIdx];
	       RDKit::Bond *bond = rdkm.getBondBetweenAtoms(idx_c, *nbrIdx);
	       if (bond) {
		  it = bond_map.find(iat);
		  if (it == bond_map.end()) {
		     // not there, add it.
		     bond_map[iat].push_back(*nbrIdx);
		  } else {
		     std::vector<int> &v = it->second;
		     if (std::find(v.begin(), v.end(), *nbrIdx) == v.end())
			v.push_back(*nbrIdx);
		  }
	       }
	    }
	    nbrIdx++;
	 }
      }
   }

   // check the bond map (debugging output)
   //
   if (0) { 
      for(it=bond_map.begin(); it!=bond_map.end(); it++) {
	 try {
	    std::string name;
	    rdkm[it->first]->getProp("name", name);
	    std::cout << "   " << it->first << " " << name <<  " : ";
	    for (unsigned int i=0; i<it->second.size(); i++)
	       std::cout << it->second[i] << " ";
	    std::cout << " :::: ";
	    for (unsigned int i=0; i<it->second.size(); i++) {
	       rdkm[it->second[i]]->getProp("name", name);
	       std::cout << name << " ";
	    }
	    std::cout << std::endl;
	 }
	 catch (const KeyErrorException &kee) {
	 }
      }
   }


   unsigned int n_max_bonds = 6; // max ring size of interest (for COD)

   // I don't think that this is needed now we have rdkit ring info.
   //
   if (false) { 
      // Return all the paths by which we get back to idx (using upto
      // (and including) n_max_bonds bonds (include backward paths also).
      // 
      for(it=bond_map.begin(); it!=bond_map.end(); it++) {
	 unsigned int idx = it->first;
	 // this returns all the paths! (in forward and backward directions)
	 std::vector<std::vector<int> > paths = trace_path(idx, bond_map, n_max_bonds);
	 if (paths.size() > 1) {
	    std::vector<int> path_sizes(paths.size());
	    for (unsigned int ii=0; ii<paths.size(); ii++)
	       path_sizes[ii] = paths[ii].size();
	    std::sort(path_sizes.begin(), path_sizes.end());

	    // dangerous interger arithmetic?
	    // we want every other path, not 5(f), 5(b), 6(f), 6(b), so
	    // sort them (smallest first) and take every other one.
	    // std::cout << "here are the path sizes: " << std::endl;
	    // for (unsigned int jj=0; jj<path_sizes.size(); jj++) { 
	    // std::cout << "   " << jj << " " << path_sizes[jj] << std::endl;
	    
	    // std::cout << "... about to do dangerous " << std::endl;
	    std::vector<int> filtered_path_sizes(path_sizes.size()/2);
	    for (unsigned int ii=0; ii<filtered_path_sizes.size(); ii++) { // every other
	       filtered_path_sizes[ii] = path_sizes[ii*2];
	       // std::cout << "set filtered_path_sizes[" << ii << "] to  path_sizes[ii*2] "
	       // << path_sizes[ii*2] << std::endl;
	    }

	    rdkm[it->first]->setProp("ring", filtered_path_sizes);
	 
	    // std::cout << "..... done dangerous" << std::endl;
	 }
      }
   }
}

bool
cod::atom_types_t::is_ring_member(unsigned int iat_ui,
				  const std::vector<std::vector<int> > &fused_rings) {

   int iat = iat_ui;
   for (unsigned int ir=0; ir<fused_rings.size(); ir++) { 
      for (unsigned int ii=0; ii<fused_rings[ir].size(); ii++) {
	 if (fused_rings[ir][ii] == iat)
	    return true;
      }
   }
   return false;
}


std::vector<std::vector<int> >
cod::atom_types_t::trace_path(unsigned int idx,
		const std::map<int, std::vector<int> > &bond_map,
		unsigned int n_max_bonds) {

   std::vector<int> in_path;
   return trace_path(idx, in_path, idx, bond_map, n_max_bonds);
}

std::vector<std::vector<int> > 
cod::atom_types_t::trace_path(unsigned int idx,
		std::vector<int> in_path_indices,
		unsigned int target_idx,
		const std::map<int, std::vector<int> > &bond_map,
		unsigned int level) {

   std::vector<std::vector<int> > vr;
   if (level == 0) {
      return vr; // empty
   } else {

      in_path_indices.push_back(idx);
      std::map<int, std::vector<int> >::const_iterator it_1 = bond_map.find(idx);
      const std::vector<int> &neighbs = it_1->second;
      for (unsigned int in=0; in<neighbs.size(); in++) {

	 bool do_recursion = false; // clumsy setting of "shall we look deeper"

	 if (neighbs[in] == int(target_idx)) {
	    if (in_path_indices.size() > 2) { // not just out and back
	       vr.push_back(in_path_indices);
	       return vr;
	    }
	 } else { 
	    if (std::find(in_path_indices.begin(), in_path_indices.end(), neighbs[in]) ==
		in_path_indices.end())
	       do_recursion = true; // only if the forward neighbor is not in current path shall
	                            // we continue deeper.
	 }
	 
	 if (do_recursion) {
	    std::vector<std::vector<int> > v = trace_path(neighbs[in],
							  in_path_indices,
							  target_idx,
							  bond_map,
							  level-1);
	    if (v.size() > 0)
	       for (unsigned int i=0; i<v.size(); i++)
		  vr.push_back(v[i]);
	 }
      }
   }
   return vr;
}



// can throw a std::runtime_error
//
// The lower the level the less
// sophisticated is the atom typing - level 0 is the lowest level atom
// type.

// nb_level is optional arg, default 0 (i.e. full atom type)
// 
// The nb_level is the neighbour level. 1 is first-neighbour
// (i.e. directly bonded) 2 is related via angle and 3 is related by
// torsion.  nb_level = 3 info is only needed for when the atom_p is
// in the ring or NB-1 of atom_p is in the ring and NB-2 is in the
// ring.
//
// We pass base because (if it is not null) we are not interested in
// tracing back up the tree (back to root via base atoms).
//
// NB-3 info should be kept separate
//
// atom_parent_p is the bonded parent of atom_p (used to check for
// back-tracking)
//
// atom_base_p: the central (nb-level 0) atom for which we are
// calculating the atom types eventually.
// 
// std::pair<std::string, std::list<cod::third_neighbour_info_t> >
cod::atom_type_t
cod::atom_types_t::get_cod_atom_type(const RDKit::Atom *atom_base_p,
				     const RDKit::Atom *atom_nb_1_p,
                                     const RDKit::Atom *atom_parent_p,
				     const RDKit::Atom *atom_p,
				     const RDKit::ROMol &rdkm,
				     int nb_level) {

   // std::cout << "get_cod_atom_type() called with nb_level " << nb_level << std::endl;

   // this could/should be in its own block
   if (nb_level == atom_type_t::COLON_DEGREE_TYPE) {
      // c.f. the atom_level_2_component_type constructor
      RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
      boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(atom_base_p);
      std::vector<unsigned int> v; // hybridizations
      while (nbrIdx != endNbrs) {
	 const RDKit::Atom *at = rdkm[*nbrIdx];
	 const RDKit::Atom *neigh_atom_p = at;
	 if (neigh_atom_p != atom_parent_p) {
	    // RDKit::Atom::HybridizationType hy = neigh_atom_p->getHybridization();
	    // int h = hybridization_to_int(hy);
	    int h = neigh_atom_p->getDegree();
	    v.push_back(h);
	 }
	 nbrIdx++;
      }
      std::sort(v.begin(), v.end());
      std::reverse(v.begin(), v.end());
      cod::atom_type_t at(v);
      return at;
   }

   atom_type_t atom_type; // returned type
      
   std::list<third_neighbour_info_t> tnil;

   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   int n = atom_p->getAtomicNum();
   std::string atom_ele = tbl->getElementSymbol(n);
   std::string atom_ring_string = "";
   // thses are data types that we can add to atom properties and are
   // linked together (by hand), i.e. whenever we add to
   // ring_size_info, we add to ring_atom_info.
   std::vector<int> ring_size_vec; // key: "ring_size"
   std::vector<int> ring_arom_vec; // key: "ring_arom"
   // std::vector<ring_info_t> ring_info_vec;

   if (nb_level != 3 && nb_level != atom_type_t::COLON_DEGREE_TYPE) {

      atom_ring_string = make_ring_info_string(atom_p).second;
      
   }

   // recursion block
   //
   if (nb_level < 3) { // this is always true, I think
      
      std::vector<atom_type_t> neighbour_types;

      RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
      boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(atom_p);
      while (nbrIdx != endNbrs) {
	 const RDKit::Atom *at = rdkm[*nbrIdx];
	 const RDKit::Atom *neigh_atom_p = at;

	 if (neigh_atom_p == atom_parent_p) {
	    // neighbour of central atom was back to parent.
	 } else {

	    // neighb_type is the type for the next level.  We don't want to
	    // add level_3 neighb_types (this function called with nb_level == 2)
	    //
	    if (nb_level < 2) {
	       // std::pair<std::string, std::list<third_neighbour_info_t> > s_pair =
	       // get_cod_atom_type(atom_base_p, atom_parent_p, atom_p, neigh_atom_p,
	       //                   rdkm, nb_level+1);
	       atom_type_t atom_type_local =
		  get_cod_atom_type(atom_base_p, atom_parent_p, atom_p, neigh_atom_p,
				    rdkm, nb_level+1);
	       
	       neighbour_types.push_back(atom_type_local);

	       if (atom_type_local.tnil.size()) {
		  std::list<third_neighbour_info_t>::const_iterator it;
		  for (it=atom_type_local.tnil.begin(); it!=atom_type_local.tnil.end(); it++)
		     tnil.push_back(*it);
		  tnil.sort();
		  tnil.unique();
	       }
	    } else {

	       third_neighbour_info_t tni =
		  get_cod_nb_3_type(atom_base_p, atom_parent_p, atom_p, neigh_atom_p, rdkm);

	       if (! tni.ele.empty()) {

		  if (false) { // debug
		     std::string name_base;
		     std::string name_this;
		     atom_base_p->getProp("name", name_base);
		     neigh_atom_p->getProp("name", name_this);
		     std::cout << "    " << name_base << " has 3rd neighb "
			       << name_this << std::endl;
		  }

		  tnil.push_back(tni);
		  tnil.sort();
		  tnil.unique();
	       }
	    }
	 }
	 nbrIdx++;
      }
      atom_ele += atom_ring_string;

      tnil.sort();
      tnil.unique();

      std::vector<std::string> nt(neighbour_types.size());
      for (unsigned int ii=0; ii<neighbour_types.size(); ii++)
	 nt[ii] = neighbour_types[ii].level_4; // FIXME

      std::pair<std::string, std::string> sp =
	 make_cod_level_3_and_4_atom_type(atom_base_p, atom_ele, nt, tnil, nb_level);

      // std::cout << "reseting atom_type " << nb_level << std::endl;
      atom_type = atom_type_t(sp.first, sp.second);
      atom_type.tnil = tnil;
   }

   if (nb_level == 0) {
      cod::atom_type_t at_for_colon = get_cod_atom_type(atom_base_p, 0, 0, 0, rdkm,
							atom_type_t::COLON_DEGREE_TYPE);
      atom_type.extract_degree_info(at_for_colon);
   }

   if (nb_level == 0) {
      // level 2 doesn't have info about the central atom, but does about the neighbours
      // 
      // std::string l2 = make_cod_level_2_atom_type(atom_base_p, rdkm);
      // atom_level_2_type l2(atom_base_p, rdkm);
      atom_type.level_2 = atom_level_2_type(atom_base_p, rdkm);
      // std::cout << "with nb_level " << nb_level << " l2 is " << l2 << std::endl;
   }

   return atom_type;
}


unsigned int 
cod::atom_types_t::get_smallest_ring_info(const RDKit::Atom *atom_p) const {

   unsigned int sr = 0;
   std::vector<int> ring_size_vec;
   try {
      atom_p->getProp("ring_size", ring_size_vec);
      if (ring_size_vec.size())
	 sr = ring_size_vec[0];
   }
   catch (const KeyErrorException &kee) {
      // no ring info on that atom, that's fine
   }
   return sr;
}


// atom_p is the current atom for which we want the 3rd-neighbour
// types, atom_base_p is the original (central, level 0) atom and we
// need that because NB-3 info is only given for atoms that share the
// ring with atom_base_p or if atom_p is a bonded neighbour of an atom
// that shares a ring with NB-3.
// 
cod::third_neighbour_info_t
cod::atom_types_t::get_cod_nb_3_type(const RDKit::Atom *atom_base_p, // the original atom 
				     const RDKit::Atom *atom_nb_1_p,
				     const RDKit::Atom *atom_parent_p,
				     const RDKit::Atom *atom_p,
				     const RDKit::ROMol &rdkm) {
   
   third_neighbour_info_t tni;

   bool in_ring = false;

   if (atom_base_p) { // should always be true

      try {
	 std::vector<int> ring_size_vec;
	 atom_base_p->getProp("ring_size", ring_size_vec);
	 in_ring = true;
      }
      catch (const KeyErrorException &kee) {
	 // no ring info on that atom, that's fine
      }
      
      if (in_ring) {

	 bool match_for_3rd_nb_info_flag =
	    check_for_3rd_nb_info(atom_base_p, atom_nb_1_p, atom_parent_p, atom_p, rdkm);
	 
	 if (match_for_3rd_nb_info_flag) {
	    
	    // i.e. in same ring or atom_p is attached to an atom in
	    // the same ring as atom_base_p

	    // what is the element and degree of atom_p?
	    
	    const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
	    int n = atom_p->getAtomicNum();
	    std::string atom_ele = tbl->getElementSymbol(n);
	    unsigned int degree = atom_p->getDegree();

	    tni = third_neighbour_info_t(atom_p, atom_ele, degree);
	 }
      }
   }
   return tni;
}

bool
cod::atom_types_t::check_for_3rd_nb_info(const RDKit::Atom *atom_base_p,
					 const RDKit::Atom *atom_nb_1_p,
					 const RDKit::Atom *atom_parent_p,
					 const RDKit::Atom *atom_p,
					 const RDKit::ROMol &rdkm) {

   // "related" means "3rd neighb info is needed for this atom"
   //
   bool related = false;

   RDKit::RingInfo* ring_info_p = rdkm.getRingInfo();
   unsigned int n_rings = ring_info_p->numRings();
   std::vector<std::vector<int> > atomRings = ring_info_p->atomRings();
   
   // are atom_parent_p and atom_base_p in the same ring?
   // 
   for (unsigned int i_ring=0; i_ring<n_rings; i_ring++) {
      const std::vector<int> &ring_atom_indices = atomRings[i_ring];
      unsigned int n_ring_atoms = ring_atom_indices.size();
      bool found_base   = false;
      bool found_parent = false;
      bool found_this   = false;
      for (unsigned int iat=0; iat<n_ring_atoms; iat++) {
	 const RDKit::Atom *ring_atom_p = rdkm[ring_atom_indices[iat]];

	 if (ring_atom_p == atom_parent_p)
	    found_parent = true;
	 if (ring_atom_p == atom_base_p)
	    found_base = true;
      }

      if (found_parent && found_base) {
	 related = true;
	 break;
      }
   }

   // or
   //
   // if NB-2 is in a ring (not necessarily the same ring as base
   // atom) then NB-3 is counted if the parent if NB-3 is also in a ring

   if (! related) {

      bool found_base   = false;
      bool found_parent = false;
      bool found_this   = false;
      bool found_nb_1   = false;
      for (unsigned int i_ring=0; i_ring<n_rings; i_ring++) {
	 const std::vector<int> &ring_atom_indices = atomRings[i_ring];
	 unsigned int n_ring_atoms = ring_atom_indices.size();
	 for (unsigned int iat=0; iat<n_ring_atoms; iat++) {
	    const RDKit::Atom *ring_atom_p = rdkm[ring_atom_indices[iat]];

	    if (ring_atom_p == atom_parent_p)
	       found_parent = true;
	    if (ring_atom_p == atom_base_p)
	       found_base = true;
	    if (ring_atom_p == atom_p)
	       found_this = true;
	    if (ring_atom_p == atom_nb_1_p)
	       found_nb_1 = true;
	 }

	 if (found_this && found_base && found_parent && found_nb_1)
	    break;
      }

      // if (found_this && found_base && found_parent && found_nb_1) {

      if (found_base && found_parent && found_nb_1) {

	 related = true;
      }
   }

   // finally exclude those atoms that are already related to base
   // atom via a (different atom set) angle.
   //
   if (related)
      if (related_via_angle(atom_base_p, atom_p, rdkm))
	 related = false;
   
   return related;
}

bool
cod::atom_types_t::related_via_angle(const RDKit::Atom *atom_in_1_p,
				     const RDKit::Atom *atom_in_2_p,
				     const RDKit::ROMol &rdkm) const {

   bool angle_related = false;
   RDKit::ROMol::ADJ_ITER nbrIdx_1, endNbrs_1;
   boost::tie(nbrIdx_1, endNbrs_1) = rdkm.getAtomNeighbors(atom_in_1_p);
   while(nbrIdx_1 != endNbrs_1) {

      const RDKit::Atom *at_mid = rdkm[*nbrIdx_1];

      RDKit::ROMol::ADJ_ITER nbrIdx_2, endNbrs_2;
      boost::tie(nbrIdx_2, endNbrs_2) = rdkm.getAtomNeighbors(at_mid);
      while(nbrIdx_2 != endNbrs_2) {

	 const RDKit::Atom *at = rdkm[*nbrIdx_2];

	 if (at == atom_in_2_p) {
	    angle_related = true;
	    break;
	 }
	 nbrIdx_2++;
      }
      if (angle_related)
	 break;
      nbrIdx_1++;
   }

   return angle_related;
}


	 


// we are making cod types for which base atom?
// base_atom_p.
//
// first is 3rd level (without 3rd neighbour info) and second is full (with neighbour info)
// 
std::pair<std::string, std::string>
cod::atom_types_t::make_cod_level_3_and_4_atom_type(const RDKit::Atom *base_atom_p,
						    const std::string &atom_ele,
						    const std::vector<std::string> &neighbour_types,
						    const std::list<third_neighbour_info_t> &tnil,
						    int nb_level) {

   std::string s;

   std::vector<std::string> n = sort_neighbours(neighbour_types, nb_level);

   if (false) { // debugging
      std::string name;
      if (base_atom_p)
	 base_atom_p->getProp("name", name);
      else
	 name = "<null-base-atom>";
      std::cout << "---- in make_cod_type() atom name: " << name
		<< "  atom_ele " << atom_ele
		<< " level " << nb_level
		<< " has sorted neighbour_types:\n";
      for (unsigned int ii=0; ii<n.size(); ii++)
	 std::cout << "   " << n[ii] << std::endl;
   
      if (tnil.size()) {
	 std::cout << " and has the following 3rd neighbours:\n";
	 std::list<third_neighbour_info_t>::const_iterator it;

	 // debug input to this function
	 std::string name_inner;
	 for (it=tnil.begin(); it!=tnil.end(); ++it) {
	    it->atom_p->getProp("name", name_inner);
	    std::cout << "__ :" << name_inner << ": "  << it->atom_p << "   \"" << it->ele << "\" "
		      << it->degree << std::endl;
	 }
      }
   } // ------------------ debug block

   if (n.size() == 0) {
      s = atom_ele;
   } else {
      s = atom_ele;
      std::vector<std::string> done_same;
      for (unsigned int i=0; i<n.size(); i++) {

	 // don't do anything if this neighbour_type is already in
	 // done_same (if it is of course, it has already been
	 // handled)
	 // 
	 std::vector<std::string>::const_iterator it =
	    std::find(done_same.begin(), done_same.end(), n[i]);
	 if (it == done_same.end()) {
	    // was not in done_same

	    std::vector<int> same; // indices of atoms with same type as this neighbour_type
	    for (unsigned int j=i+1; j<n.size(); j++)
	       if (n[j] == n[i])
		  same.push_back(j);

	    // now add neighbour info to s:
	    //
	    if (same.size() == 0) {
	       if (nb_level == 0) { 
		  s += "(";
		  s += n[i];
		  s += ")";
	       } else {
		  s += n[i];
	       }
	    } else {

	       // compression means "Fe2" instead of "FeFe" or "H3"
	       // instead of "HHH".  We don't want to do compression
	       // for "HH", because that results in "H2" - and thus
	       // doesn't make the string shorter (sigh).
	       bool do_compression = true;

	       if ((n[i].length() == 1) && (same.size() == 1))
		  do_compression = false;

	       // ----------------------- nb level 1 --------------------------------
	       // 
	       if (nb_level == 1) {
		  if (same.size() == 0) { 
		     s += n[i];
		  } else { 
		     if (do_compression) { 
			std::string save_s = s;
			s = "";
			s += save_s;
			s += n[i];
			s += coot::util::int_to_string(same.size()+1);
		     } else {
			std::string save_s = s;
			s = "";
			s += save_s;
			s += n[i];
			s += n[i];
		     }
		  }
	       }

	       // ----------------------- nb level 0 --------------------------------
	       // 
	       if (nb_level == 0) {
		  if (same.size() == 0) {
		     s += "(";
		     s += n[i];
		     s += ")";
		  } else {
		     s += "(";
		     s += n[i];
		     s+= ")"; // possibly don't do this if level is 0
		     s += coot::util::int_to_string(same.size()+1);

		  }
	       }
	    }
	    // we've done this (so that we pass over these later on in the n vector)
	    done_same.push_back(n[i]);
	 }
      }
   }

   std::string s_3rd_level = s;
   if (! s.empty()) {
      if (nb_level == 0) {
	 std::string t = make_cod_3rd_neighb_info_type(tnil);
	 s += t;
      }
   }
   return std::pair<std::string, std::string> (s_3rd_level, s);
}


// There are no repeats of atoms in tnil.
std::string
cod::atom_types_t::make_cod_3rd_neighb_info_type(const std::list<third_neighbour_info_t> &tnil) {

   // but if we have 1|C<3> and another 1|C<3> say, we do need to
   // consolidate those to 2|C<3>.

   // and we need to put C atoms before H atoms

   std::string s;

   if (tnil.size()) {

      // std::cout << "tnil size " << tnil.size() << std::endl;

      std::list<third_neighbour_info_t>::const_iterator it;
      std::map<std::string, int> type_map;
      std::map<std::string, int>::iterator it_map;
      for (it=tnil.begin(); it!=tnil.end(); it++) {
	 std::string key = it->ele;
	 key += "<";
	 key += coot::util::int_to_string(it->degree);
	 key += ">";
	 it_map = type_map.find(key);
	 if (it_map == type_map.end()) {
	    type_map[key] = 1;
	 } else {
	    it_map->second++;
	 } 
      }

      if (type_map.size() > 0) {

	 // which neighbour sorting method?
	 //
	 bool paul_sort = false;
	 bool fei_sort = true;

	 if (paul_sort) { 
	    s = "{";
	    std::map<std::string, int>::iterator it_last_value = type_map.end();
	    it_last_value--;
	    for(it_map=type_map.begin(); it_map!=type_map.end(); it_map++) {
	       s += coot::util::int_to_string(it_map->second);
	       s += "|";
	       s += it_map->first;
	       if (it_map != it_last_value)
		  s += ",";
	    }
	    s += "}";
	 }

	 if (fei_sort) {

	    std::vector<std::string> counted_neighbs; // e.g. "2|H<3>" - unsorted
	    for(it_map=type_map.begin(); it_map!=type_map.end(); it_map++) {
	       std::string sl;
	       sl += coot::util::int_to_string(it_map->second);
	       sl += "|";
	       sl += it_map->first;
	       counted_neighbs.push_back(sl);
	    }

	    if (false) {
	       std::cout << "pre-sort " << std::endl;
	       for (unsigned int jj=0; jj<counted_neighbs.size(); jj++)
		  std::cout << "   \"" << counted_neighbs[jj] << "\"";
	       std::cout << std::endl;
	    }

	    // now sort counted_neighbs
	    //
	    std::sort(counted_neighbs.begin(), counted_neighbs.end(),
		      fei_neighb_sorter);

	    if (false) {
	       std::cout << "post-sort " << std::endl;
	       for (unsigned int jj=0; jj<counted_neighbs.size(); jj++)
		  std::cout << "   " << counted_neighbs[jj];
	       std::cout << std::endl;
	    }

	    s = "{";
	    for (unsigned int ii=0; ii<counted_neighbs.size(); ii++) { 
	       s += counted_neighbs[ii];
	       if (ii < (counted_neighbs.size()-1))
		  s += ",";
	    }
	    s += "}";
	 }
      }
   }
   return s;
}

// static 
bool
cod::atom_types_t::fei_neighb_sorter(const std::string &a, const std::string &b) {

   // std::cout << "sort a: \"" << a << "\" b: \"" << b << "\"" << std::endl;
   if (a.length() > b.length()) {
      return true;
   } else {
      if (a.length() < b.length()) {
	 return false;
      } else {
	 if (a.size() == 0)
	    return true;
	 unsigned int i=0;
	 while ((i<a.length()) && (i<b.length())) {
	    // std::cout << "a: " << a << " b: " << b << " " << i << std::endl;
	    if (std::toupper(a[i]) < std::toupper(b[i]))
	       return true;
	    else
	       if (std::toupper(a[i]) >= std::toupper(b[i]))
		  return false;
	    i++;
	 }
	 return true;
      }
   }
}

// static
bool
cod::atom_types_t::neighbour_sorter(const std::string &a, const std::string &b) {
   if (a.length() > b.length())
      return true;
   if (a.length() < b.length())
      return false;
   return (a < b);
} 

std::vector<std::string>
cod::atom_types_t::sort_neighbours(const std::vector<std::string> &neighbours_in,
				   int level) {

   std::vector<std::string> n = neighbours_in;
   std::sort(n.begin(), n.end(), neighbour_sorter);
   return n;
}




// return 0,0 on failure
std::pair<unsigned int, unsigned int>
cod::atom_types_t::get_period_group(const RDKit::Atom *at) const {

   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   int n = at->getAtomicNum();
   std::string atom_ele = tbl->getElementSymbol(n);

   std::pair<unsigned int, unsigned int> period_group(0,0);
   std::map<std::string, std::pair<unsigned int, unsigned int> >::const_iterator it;
   it = element_period_group_map.find(atom_ele);
   if (it != element_period_group_map.end()) {
      period_group = it->second;
   } else {
      std::cout << "didn't find " << atom_ele << " in the map " << std::endl;
   } 
   return period_group;
}

// return 0 on failure
//
// (don't use this function if you can avoid it)
// 
unsigned int
cod::atom_types_t::make_hash_index(const RDKit::Atom *at) const {

   cod::primes primes(600000); // 12ms

   return make_hash_index(at, primes);
}


// return 0 on failure
unsigned int
cod::atom_types_t::make_hash_index(const RDKit::Atom *at, const cod::primes &primes) const {

   unsigned int hash_value = 0;

   std::pair<unsigned int, unsigned int> pg = get_period_group(at);

   unsigned int deg = at->getDegree();
   
   // what is the smallest ring of which this atom is part (0 is possible, 7 is max)?
   //
   unsigned int smallest_ring = get_smallest_ring_info(at);
   unsigned int hash_min_ring = 2;
   unsigned int ring_info = std::max(smallest_ring, hash_min_ring);
   
   bool arom = at->getIsAromatic();

   std::vector<unsigned int> pr = primes.get_primes();

   unsigned int prod = pr[arom] * pr[ring_info] * pr[8+deg] * pr[16+pg.first] * pr[24+pg.second];

   unsigned int table_length = 1000;

   hash_value = prod%table_length;

   if (false) { // debug

      const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
      int n = at->getAtomicNum();

      std::string name;
      at->getProp("name", name);
      
      // std::cout << "   " << name << " " << pg.first << " "
      // << pg.second << " hash-index " << hash_value << std::endl;
      
      std::cout << "   " << name << " ar: " << arom << " ri: " << ring_info
		<< " degree: " <<  8 + deg
		<< " per: "    << 16 + pg.first
		<< " group: "  << 24 + pg.second
		<< " prod " << prod << " hash-index " << hash_value
		<< std::endl;
   }
   return hash_value;
}


#endif // MAKE_ENHANCED_LIGAND_TOOLS
