
#include <vector>
#include "mmdb_manager.h"

#include "protein-geometry.hh"

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
      
      h_bond() {
	 donor = NULL;
	 acceptor = NULL;
	 donor_neigh = NULL;
	 acceptor_neigh = NULL;
      }
      h_bond(CAtom *d, CAtom *a) {
	 donor = d;
	 acceptor = a;
      } 
   };


   std::vector<h_bond>
   get_h_bonds(int donor_selHnd, int acceptor_selHnd, CMMDBManager *mol, protein_geometry &geom);

   // return the udd handle
   int  mark_donors_and_acceptors(int donor_selHnd, int acceptor_selHnd, CMMDBManager *mol,
				  protein_geometry &geom);


}
