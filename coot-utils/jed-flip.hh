

#ifndef JED_FLIP_HH
#define JED_FLIP_HH

#include "geometry/protein-geometry.hh"
#include "coot-coord-extras.hh"
#include "atom-tree.hh"

namespace coot {

   namespace util {

      // return a diagnostic message (set if needed)
      std::string jed_flip(int imol,
			   mmdb::Residue *residue_p, mmdb::Atom *clicked_atom,
			   bool invert_selection, 
			   protein_geometry *geom);

      std::string
      jed_flip_internal(coot::atom_tree_t &tree,
			const std::vector<dict_torsion_restraint_t> &interesting_torsions,
			const std::string &atom_name,
			int atom_idx,
			bool invert_selection);

      // return a non-null string on a problem
      //
      std::string
      jed_flip_internal(atom_tree_t &tree,
			const dict_torsion_restraint_t &torsion,
			const std::string &atom_name,
			int clicked_atom_idx,
			bool invert_selection);

   }
}

#endif // JED_FLIP_HH
