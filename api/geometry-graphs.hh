#ifndef GEOMETRY_GRAPHS_HH
#define GEOMETRY_GRAPHS_HH

#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   class geometry_graph_block_info_generic_t { 
   public:
      int imol;
      residue_spec_t residue_spec;
      atom_spec_t atom_spec;
      double distortion;
      std::string info_string;
      geometry_graph_block_info_generic_t(int imol_in, 
                                          const coot::residue_spec_t &rs,
                                          const atom_spec_t &atom_spec_in, 
                                          const double &dist_in, 
                                          const std::string &str) : imol(imol_in), residue_spec(rs), atom_spec(atom_spec_in), distortion(dist_in), info_string(str) {}
   };
}


#endif // GEOMETRY_GRAPHS_HH
