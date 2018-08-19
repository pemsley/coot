
#include <vector>
#include <mmdb2/mmdb_manager.h>

// I want a fast way of knowing if atoms are bonded
// (for non-bonded contacts)
// so instead of using a dictionary, let's try this
// by distance within the same residue and by
// hard-coded atom name comparison between residues.
//
// RNA, DNA and polypeptides

namespace coot {

   // this should only be called for atoms in the same residue
   // or atoms that are next to each other in sequence
   bool are_polymer_bonded(mmdb::Atom *at_1, mmdb::Atom *at_2);

   std::vector<std::vector<unsigned int> > make_bonds(mmdb::Manager *mol, int n_selected_atoms, int mol_atom_index_handle);

   // do I want a version that works with an atom selection?

   std::vector<std::vector<unsigned int> >
   find_1_4_connections(const std::vector<std::vector<unsigned int> > &bonds_vec);
   
}
