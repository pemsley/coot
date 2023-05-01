
#ifndef GRAPHICS_LINE_HH
#define GRAPHICS_LINE_HH

#include "Cartesian.h"

class graphics_line_t {
public:
   enum cylinder_class_t { UNK, SINGLE, DOUBLE, TRIPLE }; // so that double bonds can be drawn thinner
                                                          // than single bonds (likewise triple)
   cylinder_class_t cylinder_class;
   coot::CartesianPair positions;
   bool has_begin_cap;
   bool has_end_cap;
   // int residue_index;
   // restore this when finished
   // mmdb::Residue *residue_p; // the residue for the bond (maybe there should be 2 residues_ps? because
                             // sometimes there will be 2 residues for the same graphics_line_t.
                             // Hmm.
   int model_number; // -1 is unset
   int atom_index_1;
   int atom_index_2;
#if 0
   // default single bond constructor
   graphics_line_t(const coot::CartesianPair &p, bool b, bool e) {
      positions = p;
      has_begin_cap = b;
      has_end_cap = e;
      cylinder_class = SINGLE;
      // residue_index = -1; // unset
      //residue_p = 0;
      atom_index_1 = -1;
      atom_index_2 = -1;
      model_number = -1;
   }
#endif
   // we want atom indices now, not just the residue
   // graphics_line_t(const coot::CartesianPair &p, cylinder_class_t cc, bool b, bool e, mmdb::Residue *residue_p_in);

   graphics_line_t(const coot::CartesianPair &p, cylinder_class_t cc, bool b, bool e,
		   int model_no_in,
		   int atom_index_1_in, int atom_index_2_in) : positions(p) {
      has_begin_cap = b;
      has_end_cap = e;
      cylinder_class = cc;
      atom_index_1 = atom_index_1_in;
      atom_index_2 = atom_index_2_in;
      model_number = model_no_in;
   }
   graphics_line_t() { }
};

#endif // GRAPHICS_LINE_HH

