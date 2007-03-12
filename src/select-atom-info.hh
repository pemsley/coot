
#ifndef SELECT_ATOM_INFO_HH
#define SELECT_ATOM_INFO_HH

namespace coot { 

   class select_atom_info { 
      short int b_factor_editted;
      short int occ_editted;
   public:
      int udd;
      int molecule_number;
      std::string chain_id;
      int residue_number;
      std::string insertion_code;
      std::string atom_name;
      std::string altconf;
      float b_factor;
      float occ;
      select_atom_info(int udd_in, 
		       int molecule_number_in, 
		       const std::string &chain_id_in, 
		       int residue_number_in, 
		       const std::string &insertion_code_in, 
		       const std::string &atom_name_in, 
		       const std::string &altconf_in) { 
	udd = udd_in;
	molecule_number = molecule_number_in;
	chain_id = chain_id_in;
	residue_number = residue_number_in;
	insertion_code = insertion_code_in;
	atom_name = atom_name_in;
	altconf = altconf_in;
	//
	occ_editted = 0;
	b_factor_editted = 0; 
      }
      select_atom_info() {}
      
      // add functions that set the values of the edits
      void add_b_factor_edit(float b_factor_in) { 
	 b_factor_editted = 1; 
	 b_factor = b_factor_in; 
      }
      void add_occ_edit(float occ_in) { 
	 occ_editted = 1;
	 if ((occ_in > 1.0) || (occ_in < 0.0))
	    occ = 1.0;
	 else 
	    occ = occ_in;
      }
      short int has_b_factor_edit() const { return b_factor_editted; } 
      short int has_occ_edit() const { return occ_editted; }
      // return NULL on atom not found:
      CAtom *get_atom(CMMDBManager *mol) const {
	 CAtom *at = NULL;
	 if (mol) { 
	    int SelectionHandle = mol->NewSelection();
	    int n_atoms;
	    PCAtom *atoms;
	    mol->SelectAtoms(SelectionHandle, 0,
			     (char *) chain_id.c_str(),
			     residue_number, (char *) insertion_code.c_str(),
			     residue_number, (char *) insertion_code.c_str(),
			     "*",
			     (char *) atom_name.c_str(), "*",
			     (char *) altconf.c_str());
	    mol->GetSelIndex(SelectionHandle, atoms, n_atoms);
	    if (n_atoms > 0)
	       at = atoms[0];
	    mol->DeleteSelection(SelectionHandle);
	 }
	 return at;
      }
   };

} 

#endif // SELECT_ATOM_INFO_HH
