
#include "geometry/residue-and-atom-specs.hh"
#include "atom-tools.hh"

// So that we don't draw the atoms in the "static" molecules that have intermediate atoms.
// - should make the view more clear.
//
// maybe we should return a vector of bools, a vector that corresponds to the atoms in mol_atom_sel
//
std::set<int>
coot::atom_indices_in_other_molecule(atom_selection_container_t mol_atom_sel,
				     atom_selection_container_t moving_atom_sel) {

   std::set<int> s;
   bool debug = false;

   if (debug) {
      std::cout << "Here with moving_atom_sel SelectionHandle    "    << moving_atom_sel.SelectionHandle << std::endl;
      std::cout << "Here with moving_atom_sel UDDAtomIndexHandle "    << moving_atom_sel.UDDAtomIndexHandle << std::endl;
      std::cout << "Here with moving_atom_sel UDDOldAtomIndexHandle " << moving_atom_sel.UDDOldAtomIndexHandle << std::endl;
      
      std::cout << "Here with mol_atom_sel SelectionHandle    "    << mol_atom_sel.SelectionHandle << std::endl;
      std::cout << "Here with mol_atom_sel UDDAtomIndexHandle "    << mol_atom_sel.UDDAtomIndexHandle << std::endl;
      std::cout << "Here with mol_atom_sel UDDOldAtomIndexHandle " << mol_atom_sel.UDDOldAtomIndexHandle << std::endl;
   }
      
   for (int iat=0; iat<moving_atom_sel.n_selected_atoms; iat++) {
      mmdb::Atom *at = moving_atom_sel.atom_selection[iat];
      int idx = -1;
      if (at->GetUDData(moving_atom_sel.UDDOldAtomIndexHandle, idx) == mmdb::UDDATA_Ok) {
	 s.insert(idx);
	 if (debug)
	    std::cout << "atom_indices_in_other_molecule " << iat << " " << atom_spec_t(at) << " " << idx << std::endl;
      } else {
	 if (debug)
	    std::cout << "atom_indices_in_other_molecule " << iat << " " << atom_spec_t(at) << " idx fail" << std::endl;
      }
   }

   return s;
}
