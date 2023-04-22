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
         atom_t(const std::string &n, const std::string &e, const Cartesian &p, int fc, bool arom) :
            name(n), element(e), position(p), formal_charge(fc), aromatic(arom) {}
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

         unsigned int atom_index_1;
         unsigned int atom_index_2;
         bond_type_t bond_type;
         Cartesian centre_pos; // the position of the polygen (i.e. ring of atoms) of
                               // which this bond is a part
         bool have_centre_pos;  // was the bond from a polygon or (just) an extra-cyclic bond?
         int n_ring_atoms; // ring double bond shortening factor depends on this
         bond_t(int a1, int a2, bond_type_t bt) : atom_index_1(a1), atom_index_2(a2), bond_type(bt) {
            have_centre_pos = false;
            n_ring_atoms = 0; }
      };

      //! class molecule_t
      class molecule_t {
      public:
         bool is_valid() const { return ! (atoms.empty() || bonds.empty()); }
         std::vector<atom_t> atoms;
         std::vector<bond_t> bonds;
         void add_atom(const atom_t &a) { atoms.push_back(a); }
         void add_bond(const bond_t &b) { bonds.push_back(b); }
      };
   }

}

#endif // COOT_ATOM_HH
