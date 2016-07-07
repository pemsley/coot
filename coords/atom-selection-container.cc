
#include "atom-selection-container.hh"

mmdb::Residue *
atom_selection_container_t::get_next(mmdb::Residue *residue_in) const {

   mmdb::Residue *r = NULL;

   mmdb::Chain *chain = residue_in->GetChain();
   int this_res_no = residue_in->GetSeqNum();
   int res_no_next = this_res_no + 1;
   for (int i=0; i<n_selected_atoms; i++) { 
      if (atom_selection[i]->GetChain() == chain) {
	 // for rigor we should do some testing for insertion codes here abouts
	 // std::cout << "get_next(): comparing " << atom_selection[i]->GetSeqNum() << " "
	 // << res_no_next << std::endl;
	 if (atom_selection[i]->GetSeqNum() == res_no_next) {
	    r = atom_selection[i]->GetResidue();
	    break;
	 }
      } 
   }
   return r;
}

mmdb::Residue *
atom_selection_container_t::get_previous(mmdb::Residue *residue_in) const {

   mmdb::Residue *r = NULL;

   mmdb::Chain *chain = residue_in->GetChain();
   int this_res_no = residue_in->GetSeqNum();
   int res_no_prev = this_res_no - 1;
   for (int i=0; i<n_selected_atoms; i++) { 
      if (atom_selection[i]->GetChain() == chain) {
	 // for rigor we should do some testing for insertion codes here abouts
	 if (atom_selection[i]->GetSeqNum() == res_no_prev) {
	    r = atom_selection[i]->GetResidue();
	    break;
	 }
      } 
   }
   return r;
}
