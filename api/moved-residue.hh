
#include <vector>
#include "moved-atom.hh"

namespace coot {

   // this this the namespace for several (new) libcootapi interface classes and functions
   namespace api {

      //! a moved residue - which contains a vector of moved atoms
      class moved_residue_t {
      public:
         //! the chain-id
         std::string chain_id;
         //! the residue number
         int res_no;
         //! the insertion code
         std::string ins_code;
         //! the vector of moved atoms
         std::vector<moved_atom_t> moved_atoms;
         //! constructor
         moved_residue_t(const std::string &c, int rn, const std::string &i) : chain_id(c), res_no(rn), ins_code(i) {}
         //! add an atom to the vector of moved atoms in this residue.
         void add_atom(const moved_atom_t &mva) {moved_atoms.push_back(mva); }
      };
   }
}
