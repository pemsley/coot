
#ifndef HAVE_BONDED_PAIRS_HH
#define HAVE_BONDED_PAIRS_HH

#include <string>
#include <vector>
#include <iostream>

#include <mmdb_manager.h>

namespace coot {

   // The residues here are in order.  res_1 is comp_1 and res_2 is comp_2
   class bonded_pair_t {
   public:
      CResidue *res_1;
      CResidue *res_2;
      std::string link_type;
      bool is_fixed_first;
      bool is_fixed_second;
      bonded_pair_t(CResidue *r1, CResidue *r2, bool is_fixed_first_in, bool is_fixed_second_in,
		    const std::string &lt) {
	 res_1 = r1;
	 res_2 = r2;
	 link_type = lt;
	 is_fixed_first = is_fixed_first_in;
	 is_fixed_second = is_fixed_second_in;
      }
   };
   class bonded_pair_container_t {
   public:
      std::vector<bonded_pair_t> bonded_residues;
      bool try_add(const bonded_pair_t &bp); // check for null residues too.
      unsigned int size() const { return bonded_residues.size(); }
      bonded_pair_t operator[](unsigned int i) { return bonded_residues[i]; }
      const bonded_pair_t operator[](unsigned int i) const { return bonded_residues[i]; }
      bool linked_already_p(CResidue *r1, CResidue *r2) const;
      friend std::ostream& operator<<(std::ostream &s, bonded_pair_container_t bpc);
   };
   std::ostream& operator<<(std::ostream &s, bonded_pair_container_t bpc);

}

#endif // HAVE_BONDED_PAIRS_HH
