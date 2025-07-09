#include "geometry/residue-and-atom-specs.hh"
#include "coords/Cartesian.h"

class positioned_atom_spec_t {
    public:
    coot::atom_spec_t atom_spec;
    coot::Cartesian pos1;
    coot::Cartesian pos2;
};