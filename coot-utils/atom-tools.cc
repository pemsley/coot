
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

   if (false) {
      std::cout << "Here with moving_atom_sel SelectionHandle    "    << moving_atom_sel.SelectionHandle << std::endl;
      std::cout << "Here with moving_atom_sel UDDAtomIndexHandle "    << moving_atom_sel.UDDAtomIndexHandle << std::endl;
      std::cout << "Here with moving_atom_sel UDDOldAtomIndexHandle " << moving_atom_sel.UDDOldAtomIndexHandle << std::endl;
      
      std::cout << "Here with mol_atom_sel SelectionHandle    "    << mol_atom_sel.SelectionHandle << std::endl;
      std::cout << "Here with mol_atom_sel UDDAtomIndexHandle "    << mol_atom_sel.UDDAtomIndexHandle << std::endl;
      std::cout << "Here with mol_atom_sel UDDOldAtomIndexHandle " << mol_atom_sel.UDDOldAtomIndexHandle << std::endl;
   }

   int mol_udd_atom_index_handle = mol_atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "atom index");
   
   for (int iat=0; iat<moving_atom_sel.n_selected_atoms; iat++) {
      mmdb::Atom *at = moving_atom_sel.atom_selection[iat];
      if (debug)
	 std::cout << "debug:: in atom_indices_in_other_molecule() with moving atom " << atom_spec_t(at) << std::endl;
      int idx = -1;
      if (at->GetUDData(moving_atom_sel.UDDOldAtomIndexHandle, idx) == mmdb::UDDATA_Ok) {
	 if ((idx < mol_atom_sel.n_selected_atoms) && (idx != -1)) {
	    mmdb::Atom *at_mol = mol_atom_sel.atom_selection[idx];
	    atom_spec_t moving_atom_spec(at);
	    atom_spec_t mol_atom_spec(at_mol);
	    if (moving_atom_spec.atom_name == mol_atom_spec.atom_name) {
	       if (moving_atom_spec.res_no == mol_atom_spec.res_no) {
		  if (moving_atom_spec.chain_id == mol_atom_spec.chain_id) {
		     int udd_idx_atom;
		     if (debug)
			std::cout << "debug:: here with mol_udd_atom_index_handle " << mol_udd_atom_index_handle << std::endl;
		     if (at_mol->GetUDData(mol_udd_atom_index_handle, udd_idx_atom) == mmdb::UDDATA_Ok) {
			if (udd_idx_atom == idx) {
			   s.insert(idx);
			   if (debug)
			      std::cout << "INFO:: atom_indices_in_other_molecule " << iat << " "
					<< moving_atom_spec << " " << idx << std::endl;
			} else {
			   std::cout << "WARNING:: atom_indices_in_other_molecule() rejecting atom from set because "
				     << udd_idx_atom << " is not " << idx << std::endl;
			}
		     } else {
			std::cout << "WARNING:: atom_indices_in_other_molecule() GetUDData failure " << mol_udd_atom_index_handle
				  << std::endl;
		     }
		  } else {
		     std::cout << "WARNING:: atom_indices_in_other_molecule not same chain id " << moving_atom_spec << std::endl;
		  }
	       } else {
		  std::cout << "WARNING:: atom_indices_in_other_molecule not same res_no " << moving_atom_spec << std::endl;
	       }
	    } else {
	       std::cout << "WARNING:: atom_indices_in_other_molecule not same atom_name "
			 << moving_atom_spec << " " << mol_atom_spec << std::endl;
	    }
	 } else {
	    std::cout << "WARNING:: atom_indices_in_other_molecule - bad atom index " << idx << " "
		      << mol_atom_sel.n_selected_atoms << std::endl;
	 }
      } else {
	 if (debug)
	    std::cout << "WARNING:: atom_indices_in_other_molecule " << iat << " " << atom_spec_t(at) << " idx fail" << std::endl;
      }
   }

   return s;
}
