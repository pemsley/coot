
#include "coot-h-bonds.hh"




std::vector<coot::h_bond>
coot::get_h_bonds(int selHnd_1, int selHnd_2, CMMDBManager *mol,
		  coot::protein_geometry &geom) {

   realtype min_dist = 2.4; // H-bonds are longer than this
   realtype max_dist = 3.3; // H-bonds are shorter than this
   
   std::vector<coot::h_bond> v;

   mark_donors_and_acceptors(selHnd_1, selHnd_2, mol, geom); // using UDD data

   PPCAtom sel_1_atoms = 0;
   PPCAtom sel_2_atoms = 0;
   int n_sel_1_atoms;
   int n_sel_2_atoms;
   mat44 my_matt;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   
   mol->GetSelIndex   (selHnd_1, sel_1_atoms, n_sel_1_atoms);
   mol->GetSelIndex   (selHnd_2, sel_2_atoms, n_sel_2_atoms);

   PSContact pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;

   mol->SeekContacts(sel_1_atoms, n_sel_1_atoms,
		     sel_2_atoms, n_sel_2_atoms,
		     min_dist, max_dist,
		     0, // seqDist 0 -> also in same res.
		     pscontact, n_contacts,
		     0, &my_matt, i_contact_group);

   if (n_contacts > 0) {
      if (pscontact) {
	 for (int i_contact=0; i_contact<n_contacts; i_contact++) {
	    CAtom *at_1 = sel_1_atoms[pscontact[i_contact].id1];
	    CAtom *at_2 = sel_2_atoms[pscontact[i_contact].id2];
	 }
      }
   }

   return v;
}


// using UDD data
// 
// return the UDD handle
int
coot::mark_donors_and_acceptors(int selHnd_1, int selHnd_2, CMMDBManager *mol,
				coot::protein_geometry &geom) {

   PPCAtom sel_1_atoms = 0;
   PPCAtom sel_2_atoms = 0;
   int n_sel_1_atoms;
   int n_sel_2_atoms;
   mol->GetSelIndex   (selHnd_1, sel_1_atoms, n_sel_1_atoms);
   mol->GetSelIndex   (selHnd_2, sel_2_atoms, n_sel_2_atoms);
   int udd_h_bond_type_handle = mol->RegisterUDInteger(UDR_ATOM, "hb_type");

   for (unsigned int i=0; i<n_sel_1_atoms; i++) { 
      std::string name = sel_1_atoms[i]->name;
      std::string res_name = sel_1_atoms[i]->GetResName();
      int h_bond_type = geom.get_h_bond_type(name, res_name);
      sel_1_atoms[i]->PutUDData(udd_h_bond_type_handle, h_bond_type);
      std::cout << sel_1_atoms[i]->GetChainID() << " "
		<< sel_1_atoms[i]->GetSeqNum() << " "
		<< sel_1_atoms[i]->GetResName() << " "
		<< "name: " << name << " marked as " << h_bond_type << "\n";
   }

   if (selHnd_1 != selHnd_2) {
      for (unsigned int i=0; i<n_sel_2_atoms; i++) { 
	 std::string name = sel_2_atoms[i]->name;
	 std::string res_name = sel_2_atoms[i]->GetResName();
	 int h_bond_type = geom.get_h_bond_type(name, res_name);
	 sel_2_atoms[i]->PutUDData(udd_h_bond_type_handle, h_bond_type);
	 std::cout << sel_2_atoms[i]->GetChainID() << " "
		   << sel_2_atoms[i]->GetSeqNum() << " "
		   << sel_2_atoms[i]->GetResName() << " "
		   << "name: " << name << " marked as " << h_bond_type << "\n";
      }
   }

   return udd_h_bond_type_handle;
}
