#ifndef COOT_API_MOVED_ATOM_HH
#define COOT_API_MOVED_ATOM_HH

#include <string>

namespace coot {

   // this this the namespace for several (new) libcootapi interface classes and functions
   namespace api {

      //! a moved atom.
      //!
      //! The idea is that the user interface/front end is allowewd (of course, this is Coot) to move atoms
      //! This is the means by which the user interface send back to libcootapi the information about where
      //! the user moved the atoms.
      class moved_atom_t {
      public:
         //! the atom name
         std::string atom_name;
         //! the alt conf
         std::string alt_conf;
         //! the 3D coordinates
         float x, y, z;
         //! the index (used for fast lookup). Use -1 if unknowwn.
         int index; // for fast lookup. -1 is used for "unknown"
         //! constructor
         moved_atom_t(const std::string &a, const std::string &alt, float x_in, float y_in, float z_in) :
            atom_name(a), alt_conf(alt), x(x_in), y(y_in), z(z_in), index(-1) {}
         //! constructor
         moved_atom_t(const std::string &a, const std::string &alt, float x_in, float y_in, float z_in, int idx) :
            atom_name(a), alt_conf(alt), x(x_in), y(y_in), z(z_in), index(idx) {}
      };
   }
}

#endif // COOT_API_MOVED_ATOM_HH
