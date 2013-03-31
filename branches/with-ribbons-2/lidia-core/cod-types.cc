
#ifdef MAKE_ENTERPRISE_TOOLS

#include "coot-utils.hh"
#include "cod-types.hh"
#include "GraphMol/Rings.h"
#include "GraphMol/RingInfo.h"
#include "boost/dynamic_bitset.hpp"

// can throw a std::runtime_error
std::vector<std::string>
cod::get_cod_atom_types(RDKit::ROMol &rdkm) {

   std::vector<std::string> v;

   RDKit::RingInfo* ring_info_p = rdkm.getRingInfo();
   int n_rings = ring_info_p->numRings();

   // Maybe add a vector of ring sizes to the atoms (most atoms will
   // not have a vector added (because they are not part of rings)).
   // 
   std::vector<std::vector<int> > atomRings = ring_info_p->atomRings();
   for (unsigned int i_ring=0; i_ring<n_rings; i_ring++) {
      std::vector<int> ring_atom_indices = atomRings[i_ring];
      int n_ring_atoms = ring_atom_indices.size();

      // don't include macrocycle ring info (and the like (like cycloheptane))
      if (n_ring_atoms <= 6) { 
	 for (unsigned int iat=0; iat<n_ring_atoms; iat++) { 
	    try {
	       std::vector<int> ring_size_vec;
	       rdkm[ring_atom_indices[iat]]->getProp("ring", ring_size_vec);
	       ring_size_vec.push_back(n_ring_atoms);
	       rdkm[ring_atom_indices[iat]]->setProp("ring", ring_size_vec);
	    }
	    catch (KeyErrorException) {
	       std::vector<int> ring_size_vec(1);
	       ring_size_vec[0] = n_ring_atoms;
	       rdkm[ring_atom_indices[iat]]->setProp("ring", ring_size_vec);
	    }
	 }
      }
   }

   boost::dynamic_bitset<> fusDone(n_rings);
   int curr = 0;
   RDKit::INT_INT_VECT_MAP neighMap; // std::map<int, std::vector<int> >
   RingUtils::makeRingNeighborMap(atomRings, neighMap);
   RDKit::INT_VECT fused;
   while (curr < n_rings) {
      fused.clear();
      RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
      std::vector<std::vector<int> > fused_rings;
      fused_rings.reserve(fused.size());
      for (RDKit::INT_VECT_CI rid = fused.begin(); rid != fused.end(); ++rid) {
	 fused_rings.push_back(atomRings[*rid]);
      }

      if (fused_rings.size() > 1)
	 handle_bigger_rings_from_fused_rings(rdkm, fused_rings);
      curr++;
   }
   
   
   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm.beginAtoms(); ai!=rdkm.endAtoms(); ai++) {
      std::string s = get_cod_atom_type(0, *ai, rdkm);
      v.push_back(s);
   }
   
   if (v.size() != rdkm.getNumAtoms())
      throw std::runtime_error("mismatch size in get_cod_atom_types()");
      
   return v;
}

// ret
void
cod::handle_bigger_rings_from_fused_rings(RDKit::ROMol &rdkm,
					  const std::vector<std::vector<int> > &fused_rings) {


   // first make a connection "tree".  What is bonded to what (in the
   // rings).  bond map is a list of atom indices that are bonded to
   // the key (an atom index)
   //
   std::map<int, std::vector<int> > bond_map;
   std::map<int, std::vector<int> >::iterator it;

   for (unsigned int iat=0; iat<rdkm.getNumAtoms(); iat++) {
      if (is_ring_member(iat, fused_rings)) {
	 RDKit::Atom *this_at = rdkm[iat].get();
	 unsigned int idx_c = this_at->getIdx();
	 RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
	 boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(this_at);
	 while(nbrIdx != endNbrs) {
	    if (is_ring_member(*nbrIdx, fused_rings)) { 
	       RDKit::ATOM_SPTR at = rdkm[*nbrIdx];
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
	 catch (KeyErrorException) {
	 }
      }
   }


   unsigned int n_max_bonds = 6; // max ring size of interest (for COD)
   
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

bool
cod::is_ring_member(unsigned int iat,   const std::vector<std::vector<int> > &fused_rings) {

   for (unsigned int ir=0; ir<fused_rings.size(); ir++) { 
      for (unsigned int ii=0; ii<fused_rings[ir].size(); ii++) {
	 if (fused_rings[ir][ii] == iat)
	    return true;
      }
   }
   return false;
}

std::vector<std::vector<int> >
cod::trace_path(unsigned int idx,
		const std::map<int, std::vector<int> > &bond_map,
		unsigned int n_max_bonds) {

   std::vector<int> in_path;
   return trace_path(idx, in_path, idx, bond_map, n_max_bonds);
}

std::vector<std::vector<int> > 
cod::trace_path(unsigned int idx,
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
      std::vector<int> neighbs = it_1->second;
      for (unsigned int in=0; in<neighbs.size(); in++) {

	 bool do_recursion = false; // clumsy setting of "shall we look deeper"

	 if (neighbs[in] == target_idx) {
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
// level is an optional arg, default 2.
//
// We pass base because (if it is not null) we are not interested in
// tracing back up the tree (back to root via base atoms).
// 
std::string
cod::get_cod_atom_type(RDKit::Atom *base_atom_p,
		       RDKit::Atom *atom_p,
		       const RDKit::ROMol &rdkm,
		       int level) {

   std::string s;
   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   int n = atom_p->getAtomicNum();
   std::string atom_ele = tbl->getElementSymbol(n);
   std::string atom_ring_string = "";
   std::vector<int> ring_info;
   if (level != 0) { 
      try {
	 atom_p->getProp("ring", ring_info);
	 std::sort(ring_info.begin(), ring_info.end());
	 atom_ring_string = "[";
	 for (unsigned int i_ring=0; i_ring<ring_info.size(); i_ring++) {
	    if (i_ring > 0)
	       atom_ring_string += ",";
	    atom_ring_string += coot::util::int_to_string(ring_info[i_ring]);
	 }
	 atom_ring_string += "]";
      }
      catch (KeyErrorException kee) {
	 // no ring info on that atom, that's fine
      }
   }

   if (level > 0) {
      std::vector<std::string> neighbour_types;
      RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
      boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(atom_p);
      while (nbrIdx != endNbrs) {
	 RDKit::ATOM_SPTR at = rdkm[*nbrIdx];
	 RDKit::Atom *neigh_atom_p = at.get();

	 if (neigh_atom_p == base_atom_p) {
	    // neighbour of central atom was back to parent.
	 } else {
	    neighbour_types.push_back(get_cod_atom_type(atom_p, neigh_atom_p, rdkm, level-1));
	 }
	 nbrIdx++;
      }
      atom_ele += atom_ring_string;
      s = make_cod_type(atom_ele, neighbour_types, level);
      // std::cout << "    level " << level << " returning composite :" << s << ":" << std::endl;
   } else {
      // set s according to the element (using lowercase for second letter if needed).
      s = atom_ele;
      s += atom_ring_string;
      // std::cout << "    level " << level << " returning :" << s << ":" << std::endl;
   }
   return s;
}

bool
cod::neighbour_sorter(const std::string &a, const std::string &b) {
   if (a.length() > b.length())
      return true;
   if (a.length() < b.length())
      return false;
   if (a < b)
      return true;
   else
      return false;
} 


std::vector<std::string>
cod::sort_neighbours(const std::vector<std::string> &neighbours_in, int level) {

   std::vector<std::string> n = neighbours_in;
   std::sort(n.begin(), n.end(), neighbour_sorter);
   return n;
}


std::string
cod::make_cod_type(const std::string &atom_ele, const std::vector<std::string> &neighbour_types,
		   int level) {

   std::string s;
   std::vector<std::string> n = sort_neighbours(neighbour_types, level);

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
	       if (level == 2) { 
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

	       // ----------------------- level 1 --------------------------------
	       // 
	       if (level == 1) {
		  if (same.size() == 0) { 
		     s += n[i];
		  } else { 
		     if (do_compression) { 
			std::string save_s = s;
			s = "";
			s += save_s;
			s += n[i];
			s += coot::util::int_to_string(same.size()+1);
			s+= ""; // possibly don't do this if level is 0
		     } else {
			std::string save_s = s;
			s = "";
			s += save_s;
			s += n[i];
			s += n[i];
			s+= ""; // possibly don't do this if level is 0
		     }
		  }
	       }

	       // ----------------------- level 2 --------------------------------
	       // 
	       if (level == 2) {
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
   return s;
} 


#endif
