#ifndef SOLVENT_ACCESSIBLE_ATOM_HH
#define SOLVENT_ACCESSIBLE_ATOM_HH

#include "lidia-core/lig-build.hh"
#include "lidia-core/lbg-shared.hh"
#include <clipper/core/coords.h>

namespace pli {
   class solvent_accessible_atom_t {
   public:
      std::string atom_name;
      clipper::Coord_orth pt;
      double solvent_accessibility;
      std::vector<coot::bash_distance_t> bash_distances;
      solvent_accessible_atom_t(const std::string &at,
                                const clipper::Coord_orth &pt_in,
                                double sa) : atom_name(at), pt(pt_in), solvent_accessibility(sa) {}
      solvent_accessible_atom_t() { }
      void add_bash_dist(double d) {
         bash_distances.push_back(coot::bash_distance_t(d));
      }
      void add_unlimited() {
         coot::bash_distance_t bd;
         bash_distances.push_back(bd);
      }
   };
}


#endif // SOLVENT_ACCESSIBLE_ATOM_HH
