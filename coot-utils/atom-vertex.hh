#ifndef ATOM_VERTEX_HH
#define ATOM_VERTEX_HH

#include "mini-mol/atom-quads.hh"

namespace coot {
   class atom_vertex {

   public:
      enum connection_type_t { START, END, STANDARD, NONE };
      connection_type_t connection_type;
      std::vector<int> forward;
      std::vector<int> backward;
      std::pair<bool,atom_index_quad> torsion_quad;
      atom_vertex() {
	 connection_type = NONE;
	 torsion_quad.first = 0;
      }
   };
}


#endif // ATOM_VERTEX_HH
