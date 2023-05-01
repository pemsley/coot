#ifndef COOT_ATOM_HH
#define COOT_ATOM_HH

#include <string>
#include "coords/Cartesian.h"

namespace coot {

   namespace simple {

      //! class atom_t
      class atom_t {
      public:
         std::string name;
         std::string element;
         Cartesian position;
         int formal_charge;
         bool aromatic;
      };

      //! class bond_t
      class bond_t {
      public:

         // from lig-build.hh
         enum bond_type_t { BOND_UNDEFINED=100, SINGLE_BOND=101, DOUBLE_BOND=102,
            TRIPLE_BOND=103, AROMATIC_BOND=4, IN_BOND=104, OUT_BOND=105,
            SINGLE_OR_DOUBLE=5, SINGLE_OR_AROMATIC=6,
            DOUBLE_OR_AROMATIC=7,
            DELOC_ONE_AND_A_HALF=8,
            BOND_ANY=9 };
         
         unsigned int atom_1;
         unsigned int atom_2;
         bond_type_t bond_type;
         Cartesian centre_pos; // the position of the polygen (i.e. ring of atoms) of
                               // which this bond is a part
         bool have_centre_pos;  // was the bond from a polygon or (just) an extra-cyclic bond?
         int n_ring_atoms; // ring double bond shortening factor depends on this
      };

      //! class molecule_t
      class molecule_t {
      public:
         std::vector<atom_t> atoms;
         std::vector<atom_t> bonds;
      };

}

#endif // COOT_ATOM_HH
