
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


cod::atom_level_2_type::atom_level_2_type(RDKit::Atom *base_atom_p,
					  const RDKit::ROMol &rdkm) {

   //
   // first is int,string, number of rings and atom-ring,string
   // e.g. 2,"C[5a,6a]"  1,"N[6a]"
   // 
   // second is the  vector of hybridization of the neighbours, e.g.(2,2,3)
   // 
   std::vector<atom_level_2_component_type> components;
   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

   RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
   boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(base_atom_p);
   while(nbrIdx != endNbrs) {
      RDKit::ATOM_SPTR at_neighb = rdkm[*nbrIdx];
      unsigned int degree = at_neighb->getDegree();
      int n = at_neighb->getAtomicNum();
      std::pair<int,std::string> ring_info = make_ring_info_string(at_neighb.get());
      
      std::string atom_ele = tbl->getElementSymbol(n);

      // std::string s = atom_ele;
      // s += ring_info.second;
      // s += "-";
      
      atom_level_2_component_type c(at_neighb.get(), rdkm);
      components.push_back(c);
      
      nbrIdx++;
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

   int n = base_atom_p->getAtomicNum();

   std::string l2;

   int i_last_component = components.size() -1;
   int n_components = components.size();
   for (int i=0; i<n_components; i++) {
      l2 += components[i].element;
      l2 += components[i].ring_info_string;
      if (components[i].neighb_hybridizations.size()) {
	 l2 += "-";
	 int n_size = components[i].neighb_hybridizations.size(); // int for comparison below
	 for (int j=0; j<n_size; j++) { 
	    if (j != 0)
	       l2 += "_";
	    l2 += coot::util::int_to_string(components[i].neighb_hybridizations[j]);
	 }
	 if (i != i_last_component)
	    l2 += ":";
      }
   }

   str = l2;
}

std::ostream &
cod::operator<<(std::ostream &s,
		const atom_level_2_type::atom_level_2_component_type &c) {

   s << "{" << c.element << " " << c.number_of_rings << " "
     << c.ring_info_string << " ";
   for (unsigned int i=0; i<c.neighb_hybridizations.size(); i++)
      s << c.neighb_hybridizations[i] << " ";
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
// also maybe this is the wrong file for this function
// 
cod::atom_type_t::atom_type_t(const std::string &s1, const cod::atom_level_2_type &l2,
			      const std::string &s3, const std::string &s4) {
   level_2 = l2;
   level_3 = s3;
   level_4 = s4;

   try {
      hash_value = coot::util::string_to_int(s1);
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

   bool status = false;
   if (la.element < lb.element) {
      return true;
   } else {
      if (la.element > lb.element) {
	 return false;
      } else {

	 // same element
	 
	 bool s = (lb.neighb_hybridizations.size() <
		   la.neighb_hybridizations.size());

	 if (s) {
	    return true;
	 } else {

	    bool sr = (lb.neighb_hybridizations.size() >
		       la.neighb_hybridizations.size());

	    if (sr) {
	       return false;
	    } else {

	       // same hybridization (of 1st neighb)

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
	       
			   if (la.ring_info_string < lb.ring_info_string) {
			      return true;
			   } else {
			      if (la.ring_info_string > lb.ring_info_string) {
				 return false;
			      } else {

				 // so la.neighb_hybridizations should be the same
				 // size as lb.neighb_hybridizations.

				 unsigned int n = la.neighb_hybridizations.size();

				 for (unsigned int i=0; i<n; i++) { 
				    if (la.neighb_hybridizations[i] < lb.neighb_hybridizations[i])
				       return true;
				    if (la.neighb_hybridizations[i] > lb.neighb_hybridizations[i])
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
   return status;
}
   


std::pair<int, std::string>
cod::make_ring_info_string(RDKit::Atom *atom_p) {

   int n_rings = 0;
   std::string atom_ring_string;
   
   std::vector<int> ring_size_vec;
   std::vector<int> ring_arom_vec;

   try {
      atom_p->getProp("ring_size", ring_size_vec);
      atom_p->getProp("ring_arom", ring_arom_vec);

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
	 int arom = ring_arom_vec[i_ring];
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


cod::atom_level_2_type::atom_level_2_component_type::atom_level_2_component_type(RDKit::Atom *at, const RDKit::ROMol &rdkm) {
      
   // Version 147 atom types
   // s += coot::util::int_to_string(degree);

   // 152+ we also need info on the hybridization of the neighbours of at_neighb
   //
   std::vector<RDKit::Atom::HybridizationType> hv;
   RDKit::ROMol::ADJ_ITER nbrIdx_n, endNbrs_n;
   boost::tie(nbrIdx_n, endNbrs_n) = rdkm.getAtomNeighbors(at);
   while(nbrIdx_n != endNbrs_n) {
      RDKit::ATOM_SPTR at_neighb_n = rdkm[*nbrIdx_n];
      hv.push_back(at_neighb_n->getHybridization());
      nbrIdx_n++;
   }
   std::pair<int,std::string> ring_info = make_ring_info_string(at);

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

   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   element = tbl->getElementSymbol(at->getAtomicNum());

   neighb_hybridizations = v;

}



#endif // MAKE_ENHANCED_LIGAND_TOOLS
