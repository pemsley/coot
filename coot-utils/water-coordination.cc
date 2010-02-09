

#include "coot-coord-utils.hh"

coot::util::water_coordination_t::water_coordination_t(CMMDBManager *mol, realtype radius) {

   if (! mol)
      return; 

   realtype min_dist = 0.5;
   realtype max_dist = radius;
   
   PSContact pscontact = NULL;
   int n_contacts = 0;
   long i_contact_group = 1;
   mat44 my_matt;
   CSymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   PCAtom *water_selection = NULL;
   int n_water_atoms;
   PCAtom *atom_selection = NULL;
   int n_selected_atoms;
   int SelHnd = mol->NewSelection();
   int SelHnd_waters = mol->NewSelection();
   
   mol->SelectAtoms (SelHnd, 0, "*",
		     ANY_RES, // starting resno, an int
		     "*", // any insertion code
		     ANY_RES, // ending resno
		     "*", // ending insertion code
		     "*", // any residue name
		     "*",
		     "*", // elements
		     "*"  // alt loc.
		     );

   mol->SelectAtoms (SelHnd_waters, 0, "*",
		     ANY_RES, // starting resno, an int
		     "*", // any insertion code
		     ANY_RES, // ending resno
		     "*", // ending insertion code
		     "HOH", // residue name
		     "*",   // atname
		     "*",   // elements
		     "*"    // alt loc.
		     );
   mol->GetSelIndex(SelHnd,        water_selection, n_water_atoms);
   mol->GetSelIndex(SelHnd_waters,  atom_selection, n_selected_atoms);


   mol->SeekContacts(water_selection, n_water_atoms, 
		     atom_selection, n_selected_atoms,
		     min_dist, max_dist, // min, max distances
		     0,        // seqDist 0 -> in same res also
		     pscontact, n_contacts,
		     0, &my_matt, i_contact_group);


   // subroutine?
   if (n_contacts > 0) {
      for (int i=0; i< n_contacts; i++) {
	 
	 // e.g. the NE2 in atom_selection contacts O in an HOH in the water_selection
	 // (the test for matching altconfs id done in add_contact().
	 //
	 add_contact(water_selection[pscontact[i].id1], atom_selection[pscontact[i].id2]);
      }
   } 

   mol->DeleteSelection(SelHnd);
   mol->DeleteSelection(SelHnd_waters);
   
}


coot::util::contact_atoms_info_t::contact_atom_t::contact_atom_t(CAtom *contactor, CAtom *central_atom) {

   clipper::Coord_orth co_1(   contactor->x,    contactor->y,    contactor->z);
   clipper::Coord_orth co_2(central_atom->x, central_atom->y, central_atom->z);
   dist = clipper::Coord_orth::length(co_1, co_2);
   at = contactor;
}

// But don't add them if they are not matching alt confs
// 
void
coot::util::water_coordination_t::add_contact(CAtom *atom_contactor, CAtom *atom_central) {

   std::string alt_conf_1 = atom_contactor->altLoc;
   std::string alt_conf_2 = atom_central->altLoc;

   if ((alt_conf_1 == alt_conf_2) || (alt_conf_1 == "") || (alt_conf_2 == "")) { 

      coot::util::contact_atoms_info_t::contact_atom_t con_at(atom_contactor, atom_central);
      //
      // Now where does con_at go?
      // is atom_central already in atom_contacts?
      //
      bool found = 0;
      for (unsigned int i=0; i<atom_contacts.size(); i++) {
	 if (atom_contacts[i].matches_atom(atom_central)) {
	    atom_contacts[i].add(con_at);
	    found = 1; 
	    break;
	 } 
      }

      if (! found) {
	 coot::util::contact_atoms_info_t cai(atom_central, con_at);
	 atom_contacts.push_back(cai);
      }
   }
}


std::vector<coot::util::contact_atoms_info_t>
coot::util::water_coordination_t::get_highly_coordinated_waters(int n_contacts,  // at least n_contacts
								double dist_max) const { // within dist_max

   std::vector<coot::util::contact_atoms_info_t> v;
   // using std::vector<contact_atoms_info_t> atom_contacts;
   //
   for (unsigned int i=0; i<atom_contacts.size(); i++) {
      if (atom_contacts[i].has_contacts(n_contacts, dist_max)) {
	 v.push_back(atom_contacts[i]);
      }
   }
   return v;
}
