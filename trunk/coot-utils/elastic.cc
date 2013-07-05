

#include <iostream>
#include <map>
#include "utils/coot-utils.hh" // for random()
#include "elastic.hh"

coot::elastic_network_model_t::elastic_network_model_t(CMMDBManager *mol,
						       int SelectionHandle,
						       realtype min_dist,
						       realtype max_dist,
						       int max_n_distances) {

   if (mol) { 
      PPCAtom atom_selection = NULL;
      int n_selected_atoms;
      mol->GetSelIndex(SelectionHandle, atom_selection, n_selected_atoms);
      // std::cout << "selected " << n_selected_atoms << " atoms " << std::endl;
      
      PSContact pscontact = NULL;
      int n_contacts;
      long i_contact_group = 1;
      mat44 my_matt;
      CSymOps symm;
      for (int i=0; i<4; i++) 
	 for (int j=0; j<4; j++) 
	    my_matt[i][j] = 0.0;      
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      mol->SeekContacts(atom_selection, n_selected_atoms, 
			atom_selection, n_selected_atoms,
			min_dist, max_dist, // min, max distances
			0,        // seqDist 0 -> in same res also
			pscontact, n_contacts,
			0, &my_matt, i_contact_group);
	 
      if (pscontact) {
	 // std::cout << " Found " << n_contacts  << " contacts " << std::endl;
	 if (n_contacts > 0) { 
	    if (n_contacts <= max_n_distances) {
	       // all the contacts
	       for (unsigned int i=0; i<n_contacts; i++) { 
		  elastic_network_item_t item(atom_selection[pscontact[i].id1],
					      atom_selection[pscontact[i].id2],
					      0.1);
		  d.push_back(item);
	       }
	    } else {
	       // pick some contacts from the list of contacts.
	       //
	       // pick max_n_distances from n_contacts contacts.
	       double inv_rand_max = 1.0/double(RAND_MAX);
	       std::map<int, elastic_network_item_t> contact_indices;
	       while (contact_indices.size() < max_n_distances) {
		  int random_index = int(double(coot::util::random()) * inv_rand_max * double(n_contacts));
		  elastic_network_item_t item(atom_selection[pscontact[random_index].id1],
					      atom_selection[pscontact[random_index].id2],
					      0.1);
		  contact_indices[random_index] = item;
	       }
	       std::map<int, elastic_network_item_t>::const_iterator it;
	       for (it=contact_indices.begin(); it!=contact_indices.end(); it++)
		  d.push_back(it->second);
	    }
	 } else {
	    std::cout << "problem finding contacts" << std::endl;
	 } 
      } else {
	 std::cout << "NULL pscontact " << std::endl;
      } 
   }
   std::cout << "made " << d.size() << " distance restraints" << std::endl;
}

