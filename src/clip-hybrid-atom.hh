#ifndef CLIP_HYBRID_ATOM_HH
#define CLIP_HYBRID_ATOM_HH

#include "coords/mmdb-crystal.h"
#include <mmdb2/mmdb_manager.h>

namespace coot {
   class clip_hybrid_atom {
   public:
      mmdb::Atom *atom;
      // clipper::Coord_orth pos;
      coot::Cartesian pos;
      clip_hybrid_atom() { atom = NULL; }
      clip_hybrid_atom(mmdb::Atom *mmdb_atom_p, const coot::Cartesian &p) : atom(mmdb_atom_p), pos(p) {}
   };
}

#endif // CLIP_HYBRID_ATOM_HH
