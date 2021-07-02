#ifndef SIDE_CHAIN_HH
#define SIDE_CHAIN_HH

#include "geometry/residue-and-atom-specs.hh"
#include "geometry/protein-geometry.hh"

namespace coot {
   void do_180_degree_side_chain_flip(const residue_spec_t &spec, const std::string &alt_conf,
                                      mmdb::Manager *mol, protein_geometry *pg);

}


#endif // SIDE_CHAIN_HH
