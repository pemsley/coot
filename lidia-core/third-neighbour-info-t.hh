
#ifndef THIRD_NEIGHBOUR_INFO_T 
#define THIRD_NEIGHBOUR_INFO_T

#include "use-rdkit.hh"

namespace cod { 
   class third_neighbour_info_t {
   public:
      const RDKit::Atom *atom_p;
      std::string ele;
      unsigned int degree;
      third_neighbour_info_t() {
	 degree = 0;
      }
      third_neighbour_info_t(const RDKit::Atom *atom_p_in,
			     const std::string &e,
			     unsigned int d) {
	 atom_p = atom_p_in;
	 ele = e;
	 degree = d;
      }
      bool operator==(const third_neighbour_info_t &t) const {
	 return (t.atom_p == atom_p);
      }
      bool operator<(const third_neighbour_info_t &t) const {
	 return (t.atom_p < atom_p);
      }
   };

}

#endif // THIRD_NEIGHBOUR_INFO_T

