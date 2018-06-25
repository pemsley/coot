
// the atom_selection_container_t is not part of coot-coord-utils.hh

#include "atom-selection-container.hh"

namespace coot {

   
   // So that we don't draw the atoms in the "static" molecules that have intermediate atoms.
   // - should make the view more clear.
   //
   std::set<int> atom_indices_in_other_molecule(atom_selection_container_t mol_atom_sel,
						atom_selection_container_t moving_atom_sel);


}
