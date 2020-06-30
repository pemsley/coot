
#ifndef BONDED_QUAD_HH
#define BONDED_QUAD_HH

#include "mini-mol/atom-quads.hh"

// Bond between atoms 2 and 3, with 1 and 4 used to find the positions of the bond
//
class bonded_quad_atoms : public coot::atom_quad {
public:
   enum bond_t { NONE, SINGLE, DOUBLE, TRIPLE, DELOC };
   bonded_quad_atoms() {}
   bond_t bond_type;
};

class bonded_quad_atom_names : public coot::atom_name_quad {
public:
   enum bond_t { UNASSIGNED, SINGLE, DOUBLE, TRIPLE, DELOC };
   bonded_quad_atom_names() {}
   bonded_quad_atom_names(const std::string &atom_name_0,
                          const std::string &atom_name_1,
                          const std::string &atom_name_2,
                          const std::string &atom_name_3) : atom_name_quad(atom_name_0,
                                                                           atom_name_1,
                                                                           atom_name_2,
                                                                           atom_name_3) {
      bond_type = UNASSIGNED;
   }
   bond_t bond_type;
};


#endif // BONDED_QUAD_HH

