
#ifdef MAKE_ENTERPRISE_TOOLS

#include "coot-utils.hh"
#include "cod-types.hh"

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
   
   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm.beginAtoms(); ai!=rdkm.endAtoms(); ai++) {
      
      std::string name;
      (*ai)->getProp("name", name);
      std::cout << "\n%%%%%%%%%%%%%%%%%\n%%%%% new parent: " << name << "\n"
		<< "%%%%%%%%%%%%%%%%%" << std::endl;
      
      std::string s = get_cod_atom_type(0, *ai, rdkm);
      v.push_back(s);
   }
   if (v.size() != rdkm.getNumAtoms())
      throw std::runtime_error("mismatch size in get_cod_atom_types()");
   return v;
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

   bool debug = true;
   std::string s;
   std::vector<std::string> n = sort_neighbours(neighbour_types, level);

   if (debug) { 
      std::cout << "      --- make_cod_type called with ele: " << atom_ele << "  "
		<< neighbour_types.size() << " sorted-neighbs ";
      for (unsigned int i=0; i<neighbour_types.size(); i++)
	 std::cout << "  " << n[i];
      std::cout << "  level: " << level << std::endl;
   }
   
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
   if (debug)
      std::cout << "      --- make_cod_type returns " << s << std::endl;
      
   return s;
} 


#endif
