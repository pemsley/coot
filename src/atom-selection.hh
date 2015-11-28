
namespace coot {
   
   class atom_selection_info_t { 
   public:
      enum { UNSET, BY_STRING, BY_ATTRIBUTES }; 
      int type;
      std::string chain_id;
      int resno_start;
      int resno_end;
      std::string ins_code;
      std::string altconf;
      bool alt_conf_is_set;
      // or:
      std::string atom_selection_str;
      atom_selection_info_t(const std::string &s) { 
	 atom_selection_str = s;
	 type = BY_STRING;
	 alt_conf_is_set = 0;
      }
      atom_selection_info_t(const std::string &chain_id_in, 
			    int resno_start_in, 
			    int resno_end_in,
			    const std::string &ins_code_in) { 
	 chain_id = chain_id_in;
	 resno_start = resno_start_in;
	 resno_end = resno_end_in;
	 ins_code = ins_code_in;
	 type = BY_ATTRIBUTES;
	 alt_conf_is_set = 0;
      }
      atom_selection_info_t(const std::string &chain_id_in, 
			    int resno_start_in, 
			    int resno_end_in,
			    const std::string &ins_code_in,
			    const std::string &alt_conf_in) { 
	 chain_id = chain_id_in;
	 resno_start = resno_start_in;
	 resno_end = resno_end_in;
	 ins_code = ins_code_in;
	 type = BY_ATTRIBUTES;
	 altconf = alt_conf_in;
	 alt_conf_is_set = 1;
      }
      atom_selection_info_t() { 
	 type = UNSET;
	 alt_conf_is_set = 0;
      }

      // Return the selection handle.  It is up to the caller to
      // dispose of the atom selection, with a DeleteSelection().
      int select_atoms(mmdb::Manager *mol) const;
      void using_altconf(const std::string &altconf_in) {
	 altconf = altconf_in;
	 alt_conf_is_set = 1;
      }
      std::string name() const;
      std::string mmdb_string() const;
   };

}
