#ifndef PLAIN_ATOM_OVERLAP_HH
#define PLAIN_ATOM_OVERLAP_HH

#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   class plain_atom_overlap_t {
   public:
      int ligand_atom_index;
      atom_spec_t atom_spec_1;
      atom_spec_t atom_spec_2;
      double overlap_volume;
      double r_1;
      double r_2;
      bool is_h_bond;
      plain_atom_overlap_t() : ligand_atom_index(-1), overlap_volume(-1), r_1(-1), r_2(-1), is_h_bond(false) {}
      plain_atom_overlap_t(int ligand_atom_index_in,
                           const atom_spec_t &as_1,
                           const atom_spec_t &as_2,
                           double ov,
                           double r1,
                           double r2,
                           bool isHb) : ligand_atom_index(ligand_atom_index_in),
                                        atom_spec_1(as_1),
                                        atom_spec_2(as_2),
                                        overlap_volume(ov),
                                        r_1(r1),
                                        r_2(r2),
                                        is_h_bond(isHb) {}
   };

}

#endif // PLAIN_ATOM_OVERLAP_HH
