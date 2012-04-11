
#include <unistd.h>
#include <iostream>
#include "catcheck.hh"
#include <algorithm> // for sort

int main(int argc, char **argv) {

   if (argc > 1) { 
      std::string filename = argv[1];
      CMMDBManager *mol = get_mol(filename);
      if (mol) { 
	 water_coordination_check(mol, 3.2);
      }
   } else {
      std::cout << "Usage: " << argv[0] << " <pdb-file-name>"
		<< std::endl;
   } 
   return 0;
}

CMMDBManager *get_mol(const std::string &filename) {

   CMMDBManager *MMDBManager = new CMMDBManager();
   int err = MMDBManager->ReadCoorFile((char *)filename.c_str());
   if (err) {
      std::cout << "Error reading " << filename << std::endl;
      delete MMDBManager;
      MMDBManager = 0;
   }
   return MMDBManager;
} 


bool
contact_info_t::contacts_less(const contact_info_t &a, const contact_info_t &b) {
   return (b.contact_indices.size() < a.contact_indices.size());
} 

void water_coordination_check(CMMDBManager *mol, float max_dist) {

   // Make 2 atom selections,
   //    1) the whole molecule
   //    2) the waters
   // Then do a distance check between them:

   PCAtom *whole_protein_atom_sel;
   int SelectionHandle = mol->NewSelection();
   mol->SelectAtoms (SelectionHandle, 0, "*",
		     ANY_RES, // starting resno, an int
		     "*", // any insertion code
		     ANY_RES, // ending resno
		     "*", // ending insertion code
		     "*", // any residue name
		     "*", // atom name
		     "*", // elements
		     "*"  // alt loc.
		     );
   int nSelAtoms_all;
   mol->GetSelIndex(SelectionHandle, whole_protein_atom_sel, nSelAtoms_all);

   PCAtom *waters_atom_sel = 0;
   int SelectionHandle_waters = mol->NewSelection();
   mol->SelectAtoms (SelectionHandle_waters, 0, "*",
		     ANY_RES, // starting resno, an int
		     "*", // any insertion code
		     ANY_RES, // ending resno
		     "*", // ending insertion code
		     "WAT,HOH", // any residue name
		     "*", // atom name
		     "*", // elements
		     "*"  // alt loc.
		     );
   int nSelAtoms_waters;
   mol->GetSelIndex(SelectionHandle_waters, waters_atom_sel, nSelAtoms_waters);

   std::cout << "All Atoms: " << nSelAtoms_all << "   waters: " << nSelAtoms_waters
	     << std::endl;

   mat44 my_matt;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   long i_contact_group = 1;
   PSContact contact = NULL;
   int ncontacts;
   
   mol->SeekContacts(whole_protein_atom_sel, nSelAtoms_all,
		     waters_atom_sel, nSelAtoms_waters,
		     0.01, max_dist,
		     0,
		     contact, ncontacts, 0, &my_matt, i_contact_group);

   std::vector<contact_info_t> atom_contacts(nSelAtoms_waters);
   for (int i=0; i<atom_contacts.size(); i++) {
      atom_contacts[i].this_index = i;
   }
   
   for (int i=0; i< ncontacts; i++) {
      CAtom *at = waters_atom_sel[contact[i].id2];
      CAtom *pat = whole_protein_atom_sel[contact[i].id1];
      std::string p_ele(pat->element);
      std::string w_ele(at->element);
      if (p_ele != " C") {
	 std::cout << "Accepting contact to atom with index "
		   << contact[i].id2 << " " << at->GetAtomName()
		   << at->GetResName() << " " << at->GetSeqNum() << " "
		   << at->GetChainID() << " :" << w_ele << ": by atom "
		   << pat->GetAtomName() << " " << pat->GetResName() << " "
		   << pat->GetSeqNum() << " " << pat->GetChainID() << " :"
		   << p_ele << ": " << contact[i].id1 << std::endl;
	 atom_contacts[contact[i].id2].contact_indices.push_back(contact[i].id1);
      } else {
	 std::cout << "###Reject contact to atom with index "
		   << contact[i].id2 << " " << at->GetAtomName()
		   << at->GetResName() << " " << at->GetSeqNum() << " "
		   << at->GetChainID() << " :" << w_ele << ": by atom "
		   << pat->GetAtomName() << " " << pat->GetResName() << " "
		   << pat->GetSeqNum() << " " << pat->GetChainID() << " :"
		   << p_ele << ": " << contact[i].id1 << std::endl;
      }
   }
   std::sort(atom_contacts.begin(), atom_contacts.end(), contact_info_t::contacts_less);
   std::cout << "sorted" << std::endl;

   for (int i=0; i<atom_contacts.size(); i++) {
      std::cout << "index: " << atom_contacts[i].this_index << std::endl;
      CAtom *at = waters_atom_sel[atom_contacts[i].this_index];
      std::cout << "at: " << at << std::endl;
      std::cout << atom_contacts[i].this_index << " " << at->GetResName() << " "
		<< at->GetSeqNum() << " " << at->GetChainID() << " has "
		<< atom_contacts[i].contact_indices.size() << " contacts\n";
      for (int j=0; j<atom_contacts[i].contact_indices.size(); j++) {
	 CAtom *pat = whole_protein_atom_sel[atom_contacts[i].contact_indices[j]];
	 std::string p_ele(pat->element);
	 std::cout << "   "
		   << pat->GetAtomName() << " " << pat->GetResName() << " "
		   << pat->GetSeqNum() << " " << pat->GetChainID() << " :"
		   << p_ele << ":" << std::endl;
      }
   }
} 
