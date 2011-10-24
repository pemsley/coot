

#include "bonded-pairs.hh"

bool
coot::bonded_pair_container_t::linked_already_p(CResidue *r1, CResidue *r2) const {

   bool r = 0;
   for (unsigned int i=0; i<bonded_residues.size(); i++) {
      if (((bonded_residues[i].res_1 == r1) &&
	   (bonded_residues[i].res_2 == r2)) ||
	  ((bonded_residues[i].res_1 == r2) &&
	   (bonded_residues[i].res_2 == r1))) {
	 r = 1;
	 break;
      }
   }
   return r;
}

bool
coot::bonded_pair_container_t::try_add(const coot::bonded_pair_t &bp) {

   bool found = 0;
   for (unsigned int i=0; i<bonded_residues.size(); i++) {
      if ( (bonded_residues[i].res_1 == bp.res_1 &&
	    bonded_residues[i].res_2 == bp.res_2) ||
	   (bonded_residues[i].res_1 == bp.res_2 &&
	    bonded_residues[i].res_2 == bp.res_1) ) {
	 found = 1;
	 break;
      }
   }
   
   if (! found) {
      bonded_residues.push_back(bp);
   }
   return found;
}

std::ostream&
coot::operator<<(std::ostream &s, coot::bonded_pair_container_t bpc) {

   s << "Bonded Pair Container contains " << bpc.bonded_residues.size() << " bonded residues"
     << "\n";

   for (unsigned int i=0; i<bpc.bonded_residues.size(); i++)
      s << "   " << i << "  [:"
	<< bpc[i].link_type << ": "
	<< bpc[i].res_1->GetChainID() << " " << bpc[i].res_1->GetSeqNum() << " "
	<< bpc[i].res_1->GetInsCode() << " to " << bpc[i].res_2->GetChainID() << " "
	<< bpc[i].res_2->GetSeqNum() << " " << bpc[i].res_2->GetInsCode() << "]"
	<< "\n";

   return s; 
}

