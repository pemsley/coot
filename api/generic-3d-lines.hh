#ifndef GENERIC_3D_LINES_HH
#define GENERIC_3D_LINES_HH

#include "coords/Cartesian.h"
#include "coords/graphical-bonds-container.hh"

//! This class is for generic line set information.
//! When used for environment distances, the indexing
//! is important/useful - index 0 is for non-bonded contacts
//! and index 1 is for hydrogen bonds (or other interesting bonds).
//!
class generic_3d_lines_bonds_box_t {
public:
   //! a vector of line position pairs
   std::vector<std::vector<coot::CartesianPair> > line_segments;
   generic_3d_lines_bonds_box_t() {};
   //! the bonds gnerator makes a graphical_bonds_container and libcootapi
   //! uses that to convert to this simple container
   explicit generic_3d_lines_bonds_box_t(const graphical_bonds_container &gbc) {
      for (int icol=0; icol<gbc.num_colours; icol++) {
         std::vector<coot::CartesianPair> s;
         for (int i=0; i<gbc.bonds_[icol].num_lines; i++) {
            coot::CartesianPair egl(gbc.bonds_[icol].pair_list[i].positions);
            s.push_back(egl);
         }
         line_segments.push_back(s);
      }
   }
};


#endif // GENERIC_3D_LINES_HH
