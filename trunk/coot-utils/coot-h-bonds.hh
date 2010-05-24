
#ifndef COOT_H_BONDS_HH
#define COOT_H_BONDS_HH

#include <vector>
#include <map>
#include <algorithm>

#include "mmdb_manager.h"

#include "protein-geometry.hh"
#include "coot-coord-utils.hh"

namespace coot {

   class h_bond {
   public:
      CAtom *donor;
      CAtom *acceptor;
      CAtom *donor_neigh;
      CAtom *acceptor_neigh;
      double angle_1;  // degrees
      double angle_2;
      double dist;  // H-bond length
      bool ligand_atom_is_donor;
      
      h_bond() {
	 donor = NULL;
	 acceptor = NULL;
	 donor_neigh = NULL;
	 acceptor_neigh = NULL;
	 ligand_atom_is_donor = 1;
      }
      h_bond(CAtom *d, CAtom *a) {
	 donor = d;
	 acceptor = a;
      } 
      bool operator<(const h_bond &hb_2) const {
	 return (residue_spec_t(donor) < residue_spec_t(hb_2.donor));
      }
      bool operator==(const h_bond &hb_2) const {
	 atom_spec_t sd1(donor);
	 atom_spec_t sa1(acceptor);
	 atom_spec_t sd2(hb_2.donor);
	 atom_spec_t sa2(hb_2.acceptor);
	 return ((sd1 == sd2) && (sa1 == sa2));
      } 
   };


   class h_bonds {
      // return the udd handle
      int  mark_donors_and_acceptors(int donor_selHnd, int acceptor_selHnd, CMMDBManager *mol,
				     const protein_geometry &geom);
      
      // What is the nearest neighbour of the atoms in mol?
      // 
      std::map<CAtom *,  std::vector<std::pair<CAtom *, float> > >
      make_neighbour_map(int selHnd_1, int selHnd_2, CMMDBManager *mol);

   public:
      h_bonds() {}
      
      std::vector<h_bond>
      get(int selHnd_1, int selHnd_2, CMMDBManager *mol, const protein_geometry &geom);
      
      class atom_sorter {
	 CAtom *at;
	 coot::residue_spec_t at_res_spec;
      public:
	 atom_sorter(CAtom *at_in) { at = at_in;
	    at_res_spec = coot::residue_spec_t(at);
	 }
	 bool operator()(const std::pair<CAtom *, float> &p1,
			 const std::pair<CAtom *, float> &p2) const {
	    
	    coot::residue_spec_t n1_res_spec(p1.first);
	    coot::residue_spec_t n2_res_spec(p2.first);
	    
	    // they are both in the same residue as the donor
	    // 
	    if ((n1_res_spec == at_res_spec) && (n2_res_spec == at_res_spec)) {
	       return (p1.second < p2.second);
	    } else {
	       if (n1_res_spec == at_res_spec) { 
		  // n_1 is in (and n_2 isn't), therefore n_1 is "shorter"
		  return 1;
	       } else {
		  if (n2_res_spec == at_res_spec) { 
		     // n_2 is in (and n_1 isn't), therefore n_2 is "shorter"
		     return 0;
		  } else {
		     // neither of them are in the residue
		     return (p1.second < p2.second);
		  } 
	       } 
	    } 
	 }
      };
   };
}

#endif // COOT_H_BONDS_HH

