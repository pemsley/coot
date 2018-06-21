
#include "jed-flip.hh"
#include "atom-selection-container.hh"
#include "contact-info.hh"

// return a diagnostic message (set if needed)
std::string
coot::util::jed_flip(int imol_no, mmdb::Residue *residue_p, mmdb::Atom *atom_p,
		     bool invert_selection, coot::protein_geometry *geom) {

   std::string message;

   if (residue_p) {
      if (atom_p) {
	 std::string atom_name(atom_p->GetAtomName());
	 std::string monomer_type(residue_p->GetResName());
	 std::pair<bool, dictionary_residue_restraints_t> p =
	    geom->get_monomer_restraints(monomer_type, imol_no);
	 if (! p.first) {
	    message = "monomer type ";
	    message += monomer_type;
	    message += " not found in dictionary";
	 } else {
	    bool iht = false;  // include_hydrogen_torsions_flag
	    const std::vector<dict_torsion_restraint_t> &all_torsions =
	       p.second.get_non_const_torsions(iht);
	    if (all_torsions.size() == 0) {
	       message = "No non-const torsion for residue type ";
	       message += monomer_type;
	    } else {

	       const std::vector<std::vector<std::string> > &ring_atoms_sets = p.second.get_ligand_ring_list();

	       std::vector<dict_torsion_restraint_t> interesting_torsions;
	       for (unsigned int it=0; it<all_torsions.size(); it++) {

		  bool is_ring_torsion_flag = all_torsions[it].is_ring_torsion(ring_atoms_sets);
		  if (! all_torsions[it].is_ring_torsion(ring_atoms_sets)) {
		     if (all_torsions[it].atom_id_2_4c() == atom_name)
			interesting_torsions.push_back(all_torsions[it]);
		     if (all_torsions[it].atom_id_3_4c() == atom_name)
			interesting_torsions.push_back(all_torsions[it]);
		  }
	       }
	       if (interesting_torsions.size() == 0) {
		  message = "There are no non-CONST non-ring torsions for this atom";
	       } else {
		  mmdb::PPAtom residue_atoms = 0;
		  int n_residue_atoms;
		  residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
		  int selected_atom_idx = -1; // unset

		  for(int iat=0; iat<n_residue_atoms; iat++) {
		     mmdb::Atom *at = residue_atoms[iat];
		     if (at == atom_p) {
			selected_atom_idx = iat;
			break;
		     }
		  }

		  if (selected_atom_idx == -1) {
		     message = "Selected atom was not in the residue";
		  } else {
		     // Happy Path
		     std::string alt_conf(atom_p->altLoc);
		     atom_selection_container_t residue_asc(residue_atoms, n_residue_atoms);
		     contact_info contact = getcontacts(residue_asc, monomer_type, imol_no, geom);
		     std::vector<std::vector<int> > contact_indices =
			contact.get_contact_indices_with_reverse_contacts();

		     try {
			atom_tree_t tree(contact_indices, selected_atom_idx, residue_p, alt_conf);
			message = jed_flip_internal(tree, interesting_torsions,
						    atom_name, selected_atom_idx, invert_selection);
		     }
		     catch (const std::runtime_error &rte) {
			std::cout << "RUNTIME ERROR:: " << rte.what() << " - giving up" << std::endl;
			message = "RUNTIME ERROR:: ";
			message += rte.what();
			
		     }
		  }
	       }
	    }
	 }
      }
   }
   return message;
}




// return a non-empty string on a problem
// 
std::string
coot::util::jed_flip_internal(coot::atom_tree_t &tree,
			      const std::vector<coot::dict_torsion_restraint_t> &interesting_torsions,
			      const std::string &atom_name,
			      int atom_idx,
			      bool invert_selection) {

   std::string problem_string;
   unsigned int selected_idx = 0;
   
   if (interesting_torsions.size() > 0) {

      unsigned int best_fragment_size = 9999;
      if (interesting_torsions.size() > 1) {
	 // select the best torsion based on fragment size.
	 for (unsigned int i=0; i<interesting_torsions.size(); i++) {
	    std::string atn_1 = interesting_torsions[i].atom_id_2_4c();
	    std::string atn_2 = interesting_torsions[i].atom_id_3_4c();
	    bool reverse = false; // dummy value
	    
	    std::pair<unsigned int, unsigned int> p = tree.fragment_sizes(atn_1, atn_2, reverse);
	    if (p.first < best_fragment_size) {
	       best_fragment_size = p.first;
	       selected_idx = i;
	    } 
	    if (p.second < best_fragment_size) {
	       best_fragment_size = p.second;
	       selected_idx = i;
	    } 
	 }
      }

      problem_string = jed_flip_internal(tree, interesting_torsions[selected_idx],
					 atom_name, atom_idx, invert_selection);
   }
   return problem_string;
}

// return a non-null string on a problem
// 
std::string
coot::util::jed_flip_internal(coot::atom_tree_t &tree,
			      const coot::dict_torsion_restraint_t &torsion,
			      const std::string &atom_name,
			      int clicked_atom_idx,
			      bool invert_selection) {

   std::string problem_string;

   bool reverse = false; // reverse the moving dog<->tail fragment?

   if (invert_selection)
      reverse = true;

   std::string atn_1 = torsion.atom_id_2_4c();
   std::string atn_2 = torsion.atom_id_3_4c();

   if (torsion.atom_id_3_4c() == atom_name) {
      atn_1 = torsion.atom_id_3_4c();
      atn_2 = torsion.atom_id_2_4c();
   }

   int period = torsion.periodicity();

   if (period > 1) { 

      double angle = 360/double(period);
      std::pair<unsigned int, unsigned int> p = tree.fragment_sizes(atn_1, atn_2, false);

      if (false) {  // debug
	 std::cout << "flip this torsion: " << torsion << std::endl;
	 std::cout << "DEBUG:: jed_flip_internal() fragment sizes: " << p.first << " " << p.second
		   << std::endl;
      }

      if (p.first > p.second)
	 reverse = !reverse;

      tree.rotate_about(atn_1, atn_2, angle, reverse);
   } else {
      problem_string = "Selected torsion had a periodicity of ";
      problem_string += clipper::String(period);
   }
   return problem_string;
} 

