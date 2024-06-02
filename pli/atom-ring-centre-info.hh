#ifndef ATOM_RING_CENTRE_INFO_HH
#define ATOM_RING_CENTRE_INFO_HH

#include "lidia-core/lig-build.hh"

// trivial container for a (copy of an) atom an its ring centre (if
// it has one)
class atom_ring_centre_info_t {
public:
   lig_build::atom_t atom;
   bool has_ring_centre_flag;
   lig_build::pos_t ring_centre;
   explicit atom_ring_centre_info_t(const lig_build::atom_t &at) : atom(at) {
      has_ring_centre_flag = false;
   }
   void add_ring_centre(const lig_build::pos_t &pos) {
      ring_centre = pos;
      has_ring_centre_flag = true;
   }
};
std::ostream& operator<<(std::ostream &s, atom_ring_centre_info_t wa);



#endif // ATOM_RING_CENTRE_INFO_HH
