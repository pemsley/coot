#include <sstream>

namespace coot {
   
   class atom_spec_t {
   public:
      std::string chain;
      int resno;
      std::string insertion_code;
      std::string atom_name;
      std::string alt_conf;
      int int_user_data;
      float float_user_data;
      std::string string_user_data;
      int model_number;
      atom_spec_t() {
	 chain = "unset";
	 resno = MinInt4;
	 insertion_code = "";
	 model_number = -1;
	 int_user_data = -1;
      }
      atom_spec_t(const std::string &chain_in,
		  int resno_in,
		  const std::string &insertion_code_in,
		  const std::string &atom_name_in,
		  const std::string &alt_conf_in) {
	 chain = chain_in;
	 resno = resno_in;
	 insertion_code = insertion_code_in;
	 atom_name = atom_name_in;
	 alt_conf = alt_conf_in;
	 model_number = 1;
	 int_user_data = -1;
      }
      // This presumes at is a member of a coordinate hierarchy.
      atom_spec_t(CAtom *at) {
	 if (at) { 
	    chain          = at->GetChainID();
	    resno          = at->GetSeqNum();
	    insertion_code = at->GetInsCode();
	    model_number   = at->GetModelNum();
	    atom_name      = at->name;
	    alt_conf       = at->altLoc;
	 } else {
	    chain = "unset";
	    resno = MinInt4;
	    insertion_code = "";
	    model_number = -1;
	    int_user_data = -1;
	 }
	 int_user_data = -1; // mark as "unset" (better than not setting it)
      }
      // This presumes at is a member of a coordinate hierarchy.
      atom_spec_t(CAtom *at, const std::string &user_data_string) {
	 model_number = at->GetModelNum();
	 chain = at->GetChainID();
	 resno = at->GetSeqNum();
	 insertion_code = at->GetInsCode();
	 atom_name = at->name;
	 alt_conf = at->altLoc;
	 string_user_data = user_data_string;
      }

      bool empty() const {
	 if (resno == MinInt4)
	    return true;
	 else
	    return false;
      }

      void selectatoms(CMMDBManager *mol, int SelHnd) {
	 const char *chainid = chain.c_str();
	 const char *inscode = insertion_code.c_str();
	 const char *atname  = atom_name.c_str(); // atom name
	 const char *altconf = alt_conf.c_str();
	 
	 mol->SelectAtoms(SelHnd, 0, chainid, resno, inscode, resno, inscode,
			  "*", atname, "*", altconf);
      }

      // Presumes that atom can get to SeqNum() and InsCode()? Need
      // tested against a residue not in a hierarchy.
      bool matches_spec(CAtom *atom) const;

#ifndef THIS_IS_SWIG
      bool operator==(const atom_spec_t &matcher) const {
	 bool r = false;
	 if (matcher.model_number == model_number) { 
	    if (matcher.chain == chain) {
	       if (matcher.resno == resno) {
		  if (matcher.insertion_code == insertion_code) {
		     if (matcher.atom_name == atom_name) {
			if (matcher.alt_conf == alt_conf) {
			   r = true;
			}
		     }
		  }
	       }
	    }
	 }
	 return r;
      }
#endif 

#ifndef THIS_IS_SWIG
      // we need this if atom_spec_t are used in a std::map.
      bool operator<(const atom_spec_t &matcher) const {
	 if (matcher.empty())
	    return false; 
	 if (empty())
	    return true; 
	 if (matcher.model_number < model_number) {
	    return true;
	 } else {
	    if (matcher.chain < chain) {
	       return true;
	    } else { 
	       if (matcher.resno < resno) {
		  return true;
	       } else { 
		  if (matcher.insertion_code < insertion_code) {
		     return true;
		  } else { 
		     if (matcher.atom_name < atom_name) {
			return true;
		     } else { 
			if (matcher.alt_conf < alt_conf) {
			   return true; 
			}
		     }
		  }
	       }
	    }
	 }
	 return false;
      }
#endif // THIS_IS_SWIG      
      
#ifndef THIS_IS_SWIG      
      friend std::ostream& operator<< (std::ostream& s, const atom_spec_t &spec);
#endif // THIS_IS_SWIG      
   };
   
   bool compare_atom_specs_user_float(const atom_spec_t &a1,
				      const atom_spec_t &a2);
   bool compare_atom_specs_user_float_in_pair(const std::pair<atom_spec_t, std::string> &a,
					      const std::pair<atom_spec_t, std::string> &b);
   std::pair<atom_spec_t, atom_spec_t> link_atoms(CLink *link);

   class residue_spec_t {
   public:
      int model_number;
      std::string chain;
      int resno;
      std::string insertion_code;
      int int_user_data;
      residue_spec_t(int r) {
	 resno = r;
	 chain = "";
	 insertion_code = "";
	 int_user_data = -1;
      }
      residue_spec_t(const std::string &chain_in, int r) {
	 model_number = MinInt4;
	 resno = r;
	 chain = chain_in;
	 insertion_code = "";
	 int_user_data = -1;
      }
      residue_spec_t(int model_number_in,
		     const std::string &chain_in, int r,
		     const std::string &ins_code_in) {
	 model_number = model_number_in;
	 resno = r;
	 chain = chain_in;
	 insertion_code = ins_code_in;
	 int_user_data = -1;
      }
      residue_spec_t(const std::string &chain_in, int r,
		     const std::string &ins_code_in) {
	 model_number = MinInt4;
	 resno = r;
	 chain = chain_in;
	 insertion_code = ins_code_in;
	 int_user_data = -1;
      }
      residue_spec_t(CResidue *res) {
	 if (! res) {
	    chain = "";
	    model_number = MinInt4;
	    resno = MinInt4;
	    insertion_code = "";
	 } else { 
	    chain = res->GetChainID();
	    model_number = res->GetModelNum();
	    resno = res->GetSeqNum();
	    insertion_code = res->GetInsCode();
	 }
	 int_user_data = -1;
      } 
      residue_spec_t(const atom_spec_t &atom_spec) { 
	 model_number = MinInt4;
         chain = atom_spec.chain;
         resno = atom_spec.resno;
         insertion_code = atom_spec.insertion_code;
	 int_user_data = -1;
      }
      // This one for coot_wrap_guile
      residue_spec_t() {
	 model_number = MinInt4;
	 resno = MinInt4;
	 chain = "";
	 insertion_code = "";
	 int_user_data = -1;
      }
      bool unset_p() const {
	 bool u = true;
	 if (resno != MinInt4)
	    u = false;
	 return u;
      }
      residue_spec_t next() const {
	 residue_spec_t r = *this;
	 if (resno != MinInt4)
	    r.resno += 1;
	 return r;
      }
      residue_spec_t previous() const {
	 residue_spec_t r = *this;
	 if (resno != MinInt4)
	    r.resno -= 1;
	 return r;
      }
#ifndef THIS_IS_SWIG
      bool operator==(const residue_spec_t &matcher) const {
	 if (matcher.chain == chain) {
	    if (matcher.resno == resno) {
	       if (matcher.insertion_code == insertion_code) {
		  return 1; 
	       }
	    }
	 }
	 return 0;
      }
#endif // THIS_IS_SWIG      
      
#ifndef THIS_IS_SWIG
      bool operator<(const residue_spec_t &matcher) const{
	 if (matcher.chain == chain) {
	    if (matcher.resno == resno) {
	       if (matcher.insertion_code == insertion_code) {
		  return 0; 
	       } else {
		  if (matcher.insertion_code < insertion_code)
		     return 0;
		  else
		     return 1;
	       }
	    } else {
	       if (matcher.resno < resno)
		  return 0;
	       else
		  return 1;
	    } 
	 } else {
	    if (matcher.chain < chain)
	       return 0;
	    else
	       return 1;
	 } 
	 return 0;
      }
#endif // THIS_IS_SWIG      

      std::string format() const {
	 std::ostringstream s;
	 if (!(s << *this))
	    return "";
	 return s.str();
      }

      // return an atom selection handle for the selection in the mol
      // that matches the spec.  Caller is responsible for deleting
      // the atom selection.
      //
      // selection_key_type is typically either SKEY_NEW or SKEY_OR
      // 
      int select_atoms(CMMDBManager *mol, int selhnd, int selection_key_type);

#ifndef THIS_IS_SWIG      
      friend std::ostream& operator<< (std::ostream& s, const residue_spec_t &spec);
#endif // THIS_IS_SWIG      
   };

}
