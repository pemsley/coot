
#include "main-chain.hh"


// PDBv3 FIXME

bool
coot::is_main_chain_p(mmdb::Atom *at) { 

   std::string mol_atom_name(at->name);
   if (mol_atom_name == " N  " ||
       mol_atom_name == " C  " ||
       mol_atom_name == " CA " ||
       mol_atom_name == " H  " ||
       mol_atom_name == " HA " || // CA hydrogen
       mol_atom_name == " OXT" ||
       mol_atom_name == " O  ") {
      return true;
   } else {
      std::string res_name = at->GetResName();
      if (res_name == "GLY") {
	 if (mol_atom_name == " HA2" ||
	     mol_atom_name == " HA3") {
	    return 1;
	 }
	 
	 // Perhaps N-terminal H atom?
	 mmdb::Residue *res = at->residue;
	 if (res) {
	    if (res->isNTerminus()) {
	       if (mol_atom_name == " H1 ") return true;
	       if (mol_atom_name == " H2 ") return true;
	       if (mol_atom_name == " H3 ") return true;
	    }
	 }
      }
      return 0;
   } 
}

bool
coot::is_main_chain_or_cb_p(mmdb::Atom *at) { 

   std::string mol_atom_name(at->name);
   return is_main_chain_or_cb_p(mol_atom_name);
}

// return 0 or 1
bool
coot::is_main_chain_p(const std::string &mol_atom_name) {

   if (mol_atom_name == " N  " ||
       mol_atom_name == " C  " ||
       mol_atom_name == " H  " ||
       mol_atom_name == " CA " ||
       mol_atom_name == " CB " ||
       mol_atom_name == " HA " || // CA hydrogen
       mol_atom_name == " O  ") {
      return 1;
   } else {
      return 0;
   } 
}

// return 0 or 1
bool
coot::is_main_chain_or_cb_p(const std::string &mol_atom_name) {

   if (mol_atom_name == " N  " ||
       mol_atom_name == " C  " ||
       mol_atom_name == " H  " ||
       mol_atom_name == " CA " ||
       mol_atom_name == " OXT" ||
       mol_atom_name == " CB " ||
       mol_atom_name == " HA " || // CA hydrogen
       mol_atom_name == " O  ") {
      return 1;
   } else {
      return 0;
   } 
} 


// return 0 or 1
bool coot::is_hydrogen_p(mmdb::Atom *at) {

   std::string mol_atom_ele(at->element);
   if (mol_atom_ele == " H" ||
       mol_atom_ele == " D") {
      return 1;
   } else {
      return 0;
   }
}
